########################################################################################################################
#                                                                                                                      #
#                                                                                                                      #
#                                                                                                                      #
#                                                                                                                      #
#                                                                                                                      #
#                                                                                                                      #
########################################################################################################################

import ete3
import os
import re
import linecache
import pandas as pd
import numpy as np
import subprocess

class cd:
    """
    Context manager for changing the current working directory
    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class base_filter(object):
    '''
    Base methods to identify candidate index transfers among a pool of transfers.

    To do: add help about init params
    '''

    def __init__(self, reference_tree, gene_tree_folder, aggregate_folder, reconciliation_folder,
                 overall_tree_support_thresh=80, branch_support_thresh=95, ranger_confidence_threshold=0.9,
                 ranger_executable='/work/ranger/CorePrograms/Ranger-DTL.mac', leaves_allowed=False):

        #
        # check <reference_tree> type, if either a str or a ete3 object are being provided, are treat as such
        if type(reference_tree) is str:
            self.species_tree = ete3.Tree(reference_tree, format=1)
        else:
            self.species_tree = reference_tree.copy()

        self.overall_tree_support_thresh = overall_tree_support_thresh
        self.support_threshold           = branch_support_thresh
        self.ranger_confidence_threshold = ranger_confidence_threshold
        self.leaves_allowed              = leaves_allowed
        self.gene_tree_folder            = gene_tree_folder
        self.aggregate_folder            = aggregate_folder
        self.reconciliation_folder       = reconciliation_folder
        self.ranger                      = ranger_executable


    def match_rooting(self, reference_root, tree_to_root):
        '''
        Roots tree_to_root in the same bipartition as reference_tree. Trees must be identical besides rooting position.
        :param reference_root: tree with desired root position
        :param tree_to_root: tree whose root position should match reference_tree, RF distances between trees must be
        zero.

        :return: rooted copy of tree_to_root
        '''

        tmp_tree = tree_to_root.copy()
        #
        # traverse reference_tree root's children, should be only two, starting from the one with the smallest subtree
        for node in sorted(reference_root.children, key=len):

            #
            # if the outgroup is a terminal node (rooted in a tip)
            if node.is_leaf():
                leaf = tmp_tree.get_leaves_by_name(node.name)[0]
                tmp_tree.set_outgroup(leaf)
                return tmp_tree

            #
            # otherwise...
            else:
                #
                # check monophyly of referece_tree's outgroup within tree_to_root. If not monophyletic it means that its
                #     current root is placed within the target new root.
                is_it_monophyletic, clade_type, fucking_up = tmp_tree.check_monophyly(
                    node.get_leaf_names(),
                    'name',
                    unrooted=False
                )
                #
                # if reference_tree's outgroup is monophyletic there is nothing extra to do but to set the new outgroup
                if is_it_monophyletic:
                    equivalent = tmp_tree.get_common_ancestor(node.get_leaf_names())
                    tmp_tree.set_outgroup(equivalent)

                #
                # if not monophyletic just set the root in a random taxon not part of reference_tree's outgroup
                else:
                    tmp_tree.set_outgroup(fucking_up.pop())
                    equivalent = tmp_tree.get_common_ancestor(node.get_leaf_names())
                    tmp_tree.set_outgroup(equivalent)

                return tmp_tree

    def name_branches_as_reconciliation(self, reconciliation_file, tree):
        '''
        Name branches in gene tree with naming convention provided in RANGER-DTL reconciliation file.

        :param reconciliation_file: RANGER-DTL reconciliation file
        :param tree: gene tree used for RANGER's reconciliation

        :return:  copy of provided gene tree with named internal nodes, and list o internal nodes with duplicated names.
        '''

        tmp_tree         = tree.copy()
        branches         = re.findall('^(m\d+) = LCA\[(\S+), (\S+)\]:', reconciliation_file, re.M)
        duplicated_names = {}
        for name, leaf1, leaf2 in branches:
            node = tmp_tree.get_common_ancestor(leaf1, leaf2)
            #
            # if node is already named add it to duplicated_names list and ignore new naming.
            if node.name:
                duplicated_names[name] = node.name
                continue
            node.name = name
            node.add_feature('ranger_name', name)
        return tmp_tree, duplicated_names


    def parse_aggregated(self, group):
        '''
        Retrieve list of transfers from a RANGER-DTL aggregated reconciliation file.

        :param group:

        :return: dataframe containing all transfers, and ete3 object from the gene tree
        '''

        if not os.path.isdir('%s/%s' % (self.reconciliation_folder, group)) \
                or not os.path.isfile('%s/%s' % (self.aggregate_folder, group)):
            #
            # if there is no available reconciliation data, just return a None array
            #
            return [None, None]

        aggregated = open('%s/%s' % (self.aggregate_folder, group)).read()
        with cd('%s/%s' % (self.reconciliation_folder, group)):
            gene_tree = {'named': ete3.Tree(linecache.getline('%s.output1' % group, 8), format=1)}

        gene_tree['support'] = self.match_rooting(
            gene_tree['named'],
            ete3.Tree('%s/%s.treefile' % (self.gene_tree_folder, group))
        )
        gene_tree, duplicated_names = self.name_branches_as_reconciliation(aggregated, gene_tree['support'])
        gene_tree.add_feature('group', group)

        ufboot_distribution = [node.support for node in gene_tree.traverse() if not node.is_leaf()]
        if np.percentile(ufboot_distribution, 25) < self.overall_tree_support_thresh:
            #
            #if the overall support value of the gene tree is too low, there is no reason to assess the transfer
            #   information of the gene family
            #
            return [None, None]

        num_replicates = float(re.match('Processed (\d+) files', aggregated).group(1))

        if not self.leaves_allowed:
            transfers = re.findall('^(m\d+) = .*, Transfers = ([^0]\d+?)\], \[Most Frequent mapping --> (n\d+), \
(\d+) times\], \[Most Frequent recipient --> (n\d+), (\d+) times\].', aggregated, re.M)
        else:
            transfers = re.findall('^(m\d+) = .*, Transfers = ([^0]\d+?)\], \[Most Frequent mapping --> (\S+), \
(\d+) times\], \[Most Frequent recipient --> (\S+), (\d+) times\].', aggregated, re.M)

        selected_transfers = []
        for donor_map_name, ranger_confidence, donor_name, ranger_confidence_donor, \
            recipient_name, ranger_confidence_recipient in transfers:
            if donor_map_name in duplicated_names:
                donor_map = gene_tree.search_nodes(name=duplicated_names[donor_map_name])[0]
            else:
                donor_map = gene_tree.search_nodes(name=donor_map_name)[0]

            recipient_map_search = re.search(
                '^({children[0]}|{children[1]}).*Most Frequent mapping --> {recipient}'.format(
                    recipient=recipient_name,
                    children=[child.name for child in donor_map.children]),
                aggregated, re.M)

            if recipient_map_search:
                recipient_map_name = recipient_map_search.group(1)
                if not all([donor_name, recipient_name, donor_map_name, recipient_map_name]):
                    continue
                selected_transfers.append({
                    'donor':                       donor_name,
                    'recipient':                   recipient_name,
                    'donor_map':                   donor_map_name,
                    'recipient_map':               recipient_map_name,
                    'bipartition_support':         donor_map.support,
                    'ranger_confidence':           int(ranger_confidence) / num_replicates,
                    'ranger_confidence_donor':     int(ranger_confidence_donor) / num_replicates,
                    'ranger_confidence_recipient': int(ranger_confidence_recipient) / num_replicates
                })
        transfers_df = pd.DataFrame(selected_transfers)
        transfers_df['family'] = group
        return [transfers_df, gene_tree]


    def assess_dtl_cost(self, df):
        '''

        :param df:
        :return:
        '''
        transfer_df = df.copy()

        donor_subtrees = []
        donor_maps = []
        donor_subtree_sizes = []
        for group in transfer_df.family.unique():
            with cd('%s/%s' % (self.reconciliation_folder, group)):
                gene_tree = ete3.Tree(linecache.getline('%s.output1' % group, 8), format=1)
            for donor_map in transfer_df.loc[transfer_df.family == group,
                                             'donor_map'].unique():
                donor_maps.append([group, donor_map])
                gene_donor_branch = next(gene_tree.iter_search_nodes(name=donor_map))
                donor_subtrees.append(gene_donor_branch.write(format=9))
                donor_subtree_sizes.append(len(gene_donor_branch))

        with open('tmp_ranger.input', 'w') as out:
            out.write('%s\n' % self.species_tree.write(format=9))
            out.write('\n'.join(donor_subtrees))

        subprocess.call([
            self.ranger,
            '-q',
            '-i', 'tmp_ranger.input',
            '-o', 'tmp_ranger.output'
        ])

        for (group, donor_map), subtree_size, dtl_cost in zip(
                donor_maps,
                donor_subtree_sizes,
                re.findall('^The minimum reconciliation cost is: (\d+)',
                           open('tmp_ranger.output').read(),
                           re.M)):
            transfer_df.loc[(transfer_df.donor_map == donor_map) & (transfer_df.family == group),
                            'donor_dtl_size_ratio'] = int(dtl_cost) / subtree_size

        return transfer_df

    def name_species_tree_nodes(self, reconciliation_file):
        '''

        :param reconciliation_file:
        :return:
        '''
        species_tree_named_nodes = ete3.Tree(
            linecache.getline(
                reconciliation_file, 5
            ), format=1)

        for node in species_tree_named_nodes.traverse():
            if node.is_leaf():
                continue
            else:
                equivalent_node = self.species_tree.get_common_ancestor(node.get_leaf_names())
                if equivalent_node.get_topology_id() == node.get_topology_id():
                    equivalent_node.name = node.name
                    equivalent_node.add_feature('ranger_name', node.name)
                else:
                    print('missmatching node: %s' % node.name)


    def assess_transfer_distance(self, df):
        '''

        :param df:
        :return:
        '''
        transfer_df = df.copy()
        grouped_by_donor_recipient = transfer_df.groupby(['donor', 'recipient'])
        for donor, recipient in grouped_by_donor_recipient.groups:
            transfer_df.loc[(transfer_df.donor == donor) &
                            (transfer_df.recipient == recipient),
                            'donor_recipient_distance'] = self.species_tree.get_distance(donor, recipient)

        for donor in transfer_df.donor.unique():
            transfer_df.loc[transfer_df.donor == donor,
                            'donor_depth'] = self.species_tree.get_distance(donor, topology_only=True)
            transfer_df.loc[transfer_df.donor == donor,
                            'donor_subtree_size'] = len(
                next(self.species_tree.iter_search_nodes(ranger_name=donor))
            )

        for recipient in transfer_df.recipient.unique():
            transfer_df.loc[transfer_df.recipient == recipient,
                            'recipient_depth'] = self.species_tree.get_distance(recipient, topology_only=True)
            transfer_df.loc[transfer_df.recipient == recipient,
                            'recipient_subtree_size'] = len(
                next(self.species_tree.iter_search_nodes(ranger_name=recipient))
            )

        transfer_df['donor_depth/size_ratio'] = transfer_df['donor_depth'] / \
                                                transfer_df['donor_subtree_size']
        transfer_df['recipient_depth/size_ratio'] = transfer_df['recipient_depth'] / \
                                                    transfer_df['recipient_subtree_size']

        return transfer_df


    def map_taxonomic_level(self, df, taxa_table=None):
        '''

        :param self:
        :param df: df containing all transfers
        :param taxa_table: tab-delimited file with all taxa present in the transfer df, should be part of a \
        assembly_summmary, either genbank of refseq
        :return:
        '''
        ncbi                  = ete3.NCBITaxa()

        #
        #read taxa table and make the necessary changes
        taxa_df               = pd.read_csv(taxa_table, sep='\t')
        #
        #important: RANGER understand "_" as separation between genome and gene IDs, so we have to remove it from
        #   the taxa table in order to match data structure in the rest of the pipe
        #... ranger is also very picky about its naming convetions, tips cannot contain "."
        taxa_df['Unnamed: 0'] = taxa_df['Unnamed: 0'].apply(lambda x: x.replace('_', '').split('.')[0])
        taxa_df['accession']  = taxa_df['accession'].apply(lambda x: x.replace('_', '').split('.')[0])
        taxa_df.set_index('Unnamed: 0', inplace=True)

        #
        #this df will hold all taxonomic data for all leaves in the species tree
        taxonomy_df = pd.DataFrame()

        #
        #traverse through all leaves in the species tree and add tags regarding its taxonomic classification
        for leaf in self.species_tree.get_leaf_names():
            #
            #kelsey's naming convention wasn't the most consistent, so let's try two distinct ways to find a leaf in the
            #   taxa table
            #
            if leaf in taxa_df.index:
                node_name = taxa_df.index[taxa_df.index == leaf][0]
            elif leaf in taxa_df.accession.values:
                node_name = taxa_df.query('accession==@leaf').index[0]
            else:
                #
                #if none works, fuck it...
                continue

            #does this leaf has a valid taxid in our taxa table?
            if pd.notnull(taxa_df.loc[node_name, 'taxid']):
                #
                #sweet, what is it?
                taxid = taxa_df.loc[node_name, 'taxid']

                #
                #create a lineage dict that we can manipulate
                lineage = {taxon_rank: taxon
                           #
                           #traverse all NCBI classification available for a its taxid
                           for taxon, taxon_rank in ncbi.get_rank(ncbi.get_lineage(taxid)).items()
                           }
                #
                #add the highest level taxonomic information we have
                lineage['leaf_name'] = leaf
                #
                #I said there would be a LOT of things here...
                taxonomy_df = taxonomy_df.append(lineage, ignore_index=True)
        #
        #our most specific information about each leaf will work as our "primary key"
        taxonomy_df.set_index('leaf_name', inplace=True)

        #
        #drop all columns not related to commonly used taxonomic ranks
        to_drop = []
        for column in taxonomy_df.columns:
            if column not in ['class', 'species', 'superkingdom', 'genus',
                              'order', 'phylum', 'family', 'kingdom']:
                to_drop.append(column)
        taxonomy_df.drop(to_drop, axis='columns', inplace=True)

        transfer_df = df.copy()
        for index, row in transfer_df.iterrows():
            donor_descendants = next(
                self.species_tree.iter_search_nodes(ranger_name=row.donor)
            ).get_leaf_names()
            recipient_descendants = next(
                self.species_tree.iter_search_nodes(ranger_name=row.recipient)
            ).get_leaf_names()

            donor_taxonomy     = taxonomy_df.loc[[taxon for taxon in donor_descendants     if taxon in taxonomy_df.index]]
            recipient_taxonomy = taxonomy_df.loc[[taxon for taxon in recipient_descendants if taxon in taxonomy_df.index]]

            if not donor_taxonomy.shape[0] or not recipient_taxonomy.shape[0]:
                #
                # if no descendant has a valid taxid, ignore it...
                continue

            donor_taxonomy.dropna(    axis=1, how='any', inplace=True)
            recipient_taxonomy.dropna(axis=1, how='any', inplace=True)

            donor_taxonomy = next(donor_taxonomy.loc[:,
                                  np.invert(donor_taxonomy.T.duplicated().values)
                                  ].iterrows())[1]
            recipient_taxonomy = next(recipient_taxonomy.loc[:,
                                      np.invert(recipient_taxonomy.T.duplicated().values)
                                      ].iterrows())[1]

            common_ranks = donor_taxonomy.index.intersection(recipient_taxonomy.index)

            for rank in ['species', 'genus', 'family', 'order',
                         'class', 'phylum', 'kingdom', 'superkingdom']:
                if rank in common_ranks[donor_taxonomy[common_ranks] == recipient_taxonomy[common_ranks]]:
                    break

            transfer_df.loc[index, 'transfer_level'] = rank

        return (transfer_df)