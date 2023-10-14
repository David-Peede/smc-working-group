# Import packages.
import msprime
import numpy as np

# Define a function to simulate T_{0}.
def sim_T0(k, Ne, ploidy, seed=None):
    """
    Returns the tables, root node id, and total branch length of a tree simulated under the standard coalescent.
    
    k      -- Number of chromosomes to simulate.
    Ne     -- Effective population size.
    ploidy -- Haploid or diploid coalescent units.
    seed   -- Random seed for reporducibility.
    """
    # Simulate a tree under the standard coalescent.
    ts = msprime.sim_ancestry(
        samples=[msprime.SampleSet(k, ploidy=1)],
        population_size=Ne,
        ploidy=ploidy,
        random_seed=seed,
        discrete_genome=False,
    )
    # Make a copy of the tree-seq tables for editting.
    ts_tables = ts.dump_tables()
    # Extract the root node id.
    root_id = ts.first().root
    # Extract the total branch length.
    Lx = ts.first().total_branch_length
    return ts_tables, root_id, Lx

# Define a function to extract the current tree's information.
def extract_tree_info(node_table, edge_table, right):
    """
    Returns the edge information and the total branch length of the current tree.
    
    node_table -- Tskit node table.
    edge_table -- Tskit edge table.
    right      -- The position of the recombination event.
    """
    # Mask the edge table for the current tree.
    edge_mask = (edge_table.right == right)
    c_tree_edge_table = edge_table[edge_mask]
    # Intialize a dictionary for the current tree.
    c_tree_dicc = {}
    # Intialize the total branch length of the tree.
    Lx = 0
    # For every edge on the current tree.
    for i in range(edge_mask.sum()):
        # Intialize the subdictionary for the edge.
        c_tree_dicc[i] = {}
        # Update the pertentient information to the dictionary.
        c_tree_dicc[i]['parent'] = c_tree_edge_table.parent[i]
        c_tree_dicc[i]['child'] = c_tree_edge_table.child[i]
        c_tree_dicc[i]['upper'] = node_table.time[c_tree_dicc[i]['parent']]
        c_tree_dicc[i]['lower'] = node_table.time[c_tree_dicc[i]['child']]
        c_tree_dicc[i]['length'] = c_tree_dicc[i]['upper'] - c_tree_dicc[i]['lower']
        Lx += c_tree_dicc[i]['length']
    return c_tree_edge_table, c_tree_dicc, Lx

# Define a function to determine the distance to the next recombination event.
def draw_y(rho, Lx, ploidy):
    """
    Returns the distance to the next recombination event.
    
    rho    -- Population recombination rate.
    Lx     -- Total branch length of T_{x}.
    ploidy -- Haploid or diploid coalescent units.
    """
    # Draw y.
    y = np.random.exponential((1 / ((rho / ploidy) * Lx)))
    return y

# Define a function to determine the the lineage and age of the recombination event.
def draw_g(c_tree_dicc, Lx):
    """
    Returns the recombination event information for the current tree.
    
    c_tree_dicc      -- Dictionary of the current tree's edges.
    Lx               -- Total branch length of T_{x}.
    """
    # Compute the edge weights (ie edge_length/L(x)).
    edge_weights = [(c_tree_dicc[key]['length'] / Lx) for key in c_tree_dicc.keys()]
    # Determine which edge will have the recombinatin event.
    rec_edge_key = np.random.choice(list(c_tree_dicc.keys()), p=edge_weights)
    # Determine the age of the recombination event.
    g = np.random.uniform(c_tree_dicc[rec_edge_key]['lower'], c_tree_dicc[rec_edge_key]['upper'])
    return rec_edge_key, g

# Define a function to determine the lineage and age of the next coalescent event for the smc' model.
def draw_coal_smc_prime(node_table, c_tree_edge_table, c_tree_dicc, g, Ne, ploidy):
    """
    Returns the edge and coalescent information for the next tree.
    
    node_table        -- Tskit node table.
    c_tree_edge_table -- Tskit edge table of the current tree.
    c_tree_dicc       -- Dictionary of the current tree's edges.
    g                 -- Age of the recombination event on the current tree.
    Ne                -- Effective population size.
    ploidy            -- Haploid or diploid coalescent units.
    """
    # Determine the upper bounds of the time intervals where coalescence can occur on the current tree.
    upper_bounds = [node_table.time[i] for i in np.unique(c_tree_edge_table.parent)]
    # Sort the upper bounds by age.
    sorted_upper_bounds = sorted(upper_bounds)
    # Determine the root node id.
    root_node = np.where(node_table.time == max(sorted_upper_bounds))[0].max()
    # Intialize the lower bound of the first interval.
    c_lower_bound = 0
    # Intialize the key of the edge where the coalescent event will occur.
    coal_edge_key = None
    # For every possible coalescent interval.
    for i, c_upper_bound in enumerate(sorted_upper_bounds):
        # Determine if the recombination event occurs below the upper bound of the current interval.
        if c_upper_bound > g:
            # Determine the avaiable lineages in this interval.
            available_lineages = [
                key for key in c_tree_dicc.keys() if\
                ((c_tree_dicc[key]['upper'] >= c_upper_bound) & (c_tree_dicc[key]['lower'] <= c_lower_bound))
            ] ### YOU CAN ADD THE CONDITION (key != rec_edge_key) FOR SMC ###
            # Determine the time of the coalescent event.
            coal_time = g + np.random.exponential(((Ne * ploidy) / len(available_lineages)))
            # If the the coalescent event occurs within the current time interval.
            if c_upper_bound > coal_time > c_lower_bound:
                # Determine which edge the coalescent event occurs on.
                coal_edge_key = np.random.choice(available_lineages)
                break
            # Else, re-intialize the lower bound and move on to the next interval.
            else:
                c_lower_bound = c_upper_bound
        # Else, re-intialize the lower bound and move on to the next interval.
        else:
            c_lower_bound = c_upper_bound
    # If an edge was not found within the current tree's interval.
    if coal_edge_key == None:
        # Determine the new time of coalescences above the root.
        coal_time = node_table.time[root_node] + np.random.exponential((Ne * ploidy))
    return coal_edge_key, coal_time, root_node

# Define a function to determine the lineage and age of the next coalescent event for the smc model.
def draw_coal_smc(node_table, c_tree_edge_table, c_tree_dicc, rec_edge_key, g, Ne, ploidy):
    """
    Returns the edge and coalescent information for the next tree.
    
    node_table        -- Tskit node table.
    c_tree_edge_table -- Tskit edge table of the current tree.
    c_tree_dicc       -- Dictionary of the current tree's edges.
    rec_edge_key      -- Key of the edge with the recombination event in c_tree_dicc.
    g                 -- Age of the recombination event on the current tree.
    Ne                -- Effective population size.
    ploidy            -- Haploid or diploid coalescent units.
    """
    # Determine the upper bounds of the time intervals where coalescence can occur on the current tree.
    upper_bounds = [node_table.time[i] for i in np.unique(c_tree_edge_table.parent)]
    # Sort the upper bounds by age.
    sorted_upper_bounds = sorted(upper_bounds)
    # Determine the root node id.
    root_node = np.where(node_table.time == max(sorted_upper_bounds))[0].max()
    # Intialize the lower bound of the first interval.
    c_lower_bound = 0
    # Intialize the key of the edge where the coalescent event will occur.
    coal_edge_key = None
    # For every possible coalescent interval.
    for i, c_upper_bound in enumerate(sorted_upper_bounds):
        # Determine if the recombination event occurs below the upper bound of the current interval.
        if c_upper_bound > g:
            # Determine the avaiable lineages in this interval.
            available_lineages = [
                key for key in c_tree_dicc.keys() if\
                ((c_tree_dicc[key]['upper'] >= c_upper_bound)\
                & (c_tree_dicc[key]['lower'] <= c_lower_bound)\
                & (key != rec_edge_key))
            ]
            # If there are avaiable lineages.
            if len(available_lineages) > 0:
                # Determine the time of the coalescent event.
                coal_time = g + np.random.exponential(((Ne * ploidy) / len(available_lineages)))
            # Else set the coalescent event to a variable that will fail.
            else:
                coal_time = -1
            # If the the coalescent event occurs within the current time interval.
            if c_upper_bound > coal_time > c_lower_bound:
                # Determine which edge the coalescent event occurs on.
                coal_edge_key = np.random.choice(available_lineages)
                break
            # Else, re-intialize the lower bound and move on to the next interval.
            else:
                c_lower_bound = c_upper_bound
        # Else, re-intialize the lower bound and move on to the next interval.
        else:
            c_lower_bound = c_upper_bound
    # If an edge was not found within the current tree's interval.
    if coal_edge_key == None:
        # Determine the new time of coalescences above the root.
        coal_time = node_table.time[root_node] + np.random.exponential((Ne * ploidy))
    return coal_edge_key, coal_time, root_node

# Define a function to perform a subtree pruning and regrafting that is compatible with tskit.
def ts_spr(
    k,
    node_table,
    edge_table,
    c_tree_edge_table,
    c_tree_dicc,
    rec_edge_key,
    coal_edge_key,
    coal_time,
    right,
    root_node,
    node_id,
):
    """
    Returns the updated node and edge table after performing the SPR algorithim.
    
    k                 -- Number of sampled chromosomes.
    node_table        -- Tskit node table.
    edge_table        -- Tskit edge table.
    c_tree_edge_table -- Tskit edge table of the current tree.
    c_tree_dicc       -- Dictionary of the current tree's edges.
    rec_edge_key      -- Key of the edge with the recombination event in c_tree_dicc.
    coal_edge_key     -- Key of the edge with the coalescent event in c_tree_dicc.
    coal_time         -- Time of the next coalescent event.
    right             -- The position of the recombination event.
    root_node         -- ID of the root node on the current tree.
    node_id           -- ID of the next node for Tskit compatability.
    """
    ### HIDDEN RECOMBINATION SCENARIO ###
    # If the coalescent event is hidden (ie recombination and coalesence occur on the same branch).
    if rec_edge_key == coal_edge_key:
        # Intialize a variable to indicate there was a hidden recombination event.
        hidden = 1
        # For every edge on the current tree.
        for key in c_tree_dicc.keys():
            # Extract the row index.
            row_idx = np.where(
                (edge_table.parent == c_tree_dicc[key]['parent']) & (edge_table.child == c_tree_dicc[key]['child'])
            )[0][0]
            # Update the edge table.
            edge_table[row_idx] = edge_table[row_idx].replace(right=1.0)
    # Else, the coalescent event is not hidden.
    else:
        # Intialize a variable to indicate there was no hidden recombination event.
        hidden = 0
        ### INTIALIZING TREE INFORMATION FOR SPR ###
        # Intialize a dictionary of parent nodes and children for the current tree.
        c_parent_dicc = {}
        # For every unique parent.
        for parent in np.unique(c_tree_edge_table.parent):
            # Intialize the entry
            c_parent_dicc[parent] = {
                'time': node_table.time[parent],
                'children': c_tree_edge_table[c_tree_edge_table.parent == parent].child,
            }
        # Idenitfy the broken node (ie the parent node directly above g).
        broken_node = c_tree_dicc[rec_edge_key]['parent']
        # Identify the node to be inherited (ie the parent node directly below g).
        inherited_node = c_tree_dicc[rec_edge_key]['child']
        # Idenitfy the lonely node not inherited (ie the child node of the broken node not inherited).
        lonely_node = np.setdiff1d(c_parent_dicc[broken_node]['children'], inherited_node)[0]
        # If the new coalescent event is above the current root node.
        if coal_edge_key == None:
            # Set the root node as the node below the new coalescent event.
            below_node = root_node
        # Else, the new coalescent event occurs on a current edge.
        else:
            # Identify the node below the new coalescent event.
            below_node = c_tree_dicc[coal_edge_key]['child']
        ### SPR ###
        # Intialize a dictionary to store the next tree's information.
        n_parent_dicc = c_parent_dicc.copy()
        # If the broken node is the root node and the coalescent event is above the root.
        if (broken_node == root_node) & (coal_edge_key == None):
            # Set the new coalescent as the parent node with the inherited and lonely nodes as its chidren.
            n_parent_dicc[None] = {
                'time': coal_time,
                'children': np.array([inherited_node, lonely_node], dtype='int32'),
            }
        # Else if, the broken node is the node bewlow the coalescent event.
        elif (broken_node == below_node) & (coal_edge_key != None):
            # Set the new coalescent as the parent node with the inherited and lonely nodes as its chidren.
            n_parent_dicc[None] = {
                'time': coal_time,
                'children': np.array([inherited_node, lonely_node], dtype='int32'),
            }
        # Else.
        else:
            # Set the new coalescent as the parent node with the inherited and below nodes as its chidren.
            n_parent_dicc[None] = {
                'time': coal_time,
                'children': np.array([inherited_node, below_node], dtype='int32'),
            }
        # If the below node has a parent node on the current tree.
        if (below_node != root_node) & (below_node != lonely_node):
            # Replace the below node with the new coalescent event node in the children's set.
            n_parent_dicc[c_tree_dicc[coal_edge_key]['parent']]['children'] = np.where(
                c_parent_dicc[c_tree_dicc[coal_edge_key]['parent']]['children'] == below_node,
                None, c_parent_dicc[c_tree_dicc[coal_edge_key]['parent']]['children'],
            )
        # If the broken node has a parent node on the current tree.
        if (broken_node != root_node) & (below_node != lonely_node):
            # Determine the parent node of broken node.
            broken_node_parent = c_tree_edge_table[np.where((c_tree_edge_table.child == broken_node))[0][0]].parent
            # Replace the broken node as the lonely node in the parent's node children set.
            n_parent_dicc[broken_node_parent]['children'] = np.where(
                c_parent_dicc[broken_node_parent]['children'] == broken_node,
                lonely_node, c_parent_dicc[broken_node_parent]['children'],
            )
        # Else, if broken node is not the root node.
        elif (broken_node != root_node):
            # Determine the parent node of broken node.
            broken_node_parent = c_tree_edge_table[np.where((c_tree_edge_table.child == broken_node))[0][0]].parent
            # Replace the broken node as the lonely node in the parent's node children set.
            n_parent_dicc[broken_node_parent]['children'] = np.where(
                c_parent_dicc[broken_node_parent]['children'] == broken_node,
                None, c_parent_dicc[broken_node_parent]['children'],
            )
        # Remove the broken node from the next tree.
        del n_parent_dicc[broken_node]
        ### RE-FORMATTING NODE IDS FOR TSKIT ###
        # Intialize a node id dictionary.
        node_id_dicc = {}
        # Add the node id to the node table.
        node_table.add_row(time=n_parent_dicc[None]['time'], population=0)
        # Update the node id dictionary.
        node_id_dicc[None] = node_id
        # Move the node_id counter forward.
        node_id += 1
        # If the broken node is the parent of the lonely node and the lonely node is not a leaf.
        if (c_tree_edge_table[c_tree_edge_table.child == lonely_node][0].parent == broken_node) & (lonely_node >= k):
            # Update the lonely node id so it is tskit compatible.
            node_table.add_row(time=node_table.time[lonely_node], population=0)
            # Update the node id dictionary.
            node_id_dicc[lonely_node] = node_id
            # Move the node_id counter forward.
            node_id += 1
        # Else-if the below node is the root and the lonely node is a leaf.
        elif (below_node == root_node) & (lonely_node < k):
            # Update the root node id so it is tskit compatible.
            node_table.add_row(time=node_table.time[root_node], population=0)
            # Update the node id dictionary.
            node_id_dicc[root_node] = node_id
            # Move the node_id counter forward.
            node_id += 1
        # For every parent node in the new tree.
        for key in list(n_parent_dicc.keys()):
            # If the node id needs to be updated.
            if key in node_id_dicc:
                # Update the parent node id.
                n_parent_dicc[node_id_dicc[key]] = n_parent_dicc.pop(key)
        # For every new parent node in the new tree.
        for key in n_parent_dicc.keys():
            # Update the children node ids if necessary.
            n_parent_dicc[key]['children'] = np.array(
                [node_id_dicc.get(child, child) for child in n_parent_dicc[key]['children']], dtype='int32',
            )
        ### UPDATING TSKIT TABLES ###
        # Sort the remaining nodes on the new tree by their age.
        sorted_parent_nodes = sorted(n_parent_dicc.keys(), key=lambda x: n_parent_dicc[x]['time'], reverse=True)
        # Intialize a list of new node ids.
        new_node_ids = list(node_id_dicc.values())
        # For the remaining new parent nodes.
        for n_parent in sorted_parent_nodes:
            # If the parent node is a new id.
            if n_parent in new_node_ids:
                # For every child.
                for child in n_parent_dicc[n_parent]['children']:
                    # Add the new edge to the edge table.
                    edge_table.add_row(
                        left=right, right=1.0,
                        parent=n_parent, child=child,
                    )
            # Else, the parent node is an existing node id.
            else:
                # For every child.
                for child in n_parent_dicc[n_parent]['children']:
                    # If the child node has a new node id.
                    if child in new_node_ids:
                        # Add the new edge to the edge table.
                        edge_table.add_row(
                            left=right, right=1.0,
                            parent=n_parent, child=child,
                        )
                    # Else, the child node id is an existing node.
                    else:
                        # Extract the row index from the complete edge table.
                        row_idx = np.where((edge_table.parent == n_parent) & (edge_table.child == child))[0]
                        # If, this is an existing edge on the current tree.
                        if row_idx.size > 0:
                            # Update the edge table.
                            edge_table[row_idx[0]] = edge_table[row_idx[0]].replace(right=1.0)
                        # Else, this edge is not on the current tree.
                        else:
                            # Add the new edge to the edge table.
                            edge_table.add_row(
                                left=right, right=1.0,
                                parent=n_parent, child=child,
                            )
    return node_table, edge_table, node_id, hidden

# Define a function to simulate a tree-sequence under the SMC model.
def sim_ts_smc(k, Ne, rho, ploidy, seed=None, export=False, path='./smc_ts'):
    """
    Returns a tree-sequence simulated using the Sequentially Markovian Coalescent (SMC).
    
    k      -- Number of chromosomes to simulate.
    Ne     -- Effective population size.
    rho    -- Population recombination rate.
    ploidy -- Haploid or diploid coalescent units.
    seed   -- Random seed for simulating T_{0}.
    export -- Do you want to export the tree-sequence?
    path   -- Path to export the tree-sequence. 
    """
    ## (1) Intialize the first tree, T(x)=T_{0}, at position x=0, and compute the total branch length L(x)=L_{0}. ##

    # Intialize the start position.
    x = 0
    # Simulate a tree (T_{0}) under the standard coalescent at point x=0.
    smc_tables, root_id, L0 = sim_T0(k=k, Ne=Ne, ploidy=ploidy, seed=seed)
    # Extract the tables that we will need to edit.
    node_table = smc_tables.nodes
    edge_table = smc_tables.edges
    # Intialize the next node id for tree-seq compatability.
    node_id = root_id + 1

    ## (2) Generate the distance, y=exp[(rho/2)L(x)], to the next recombination event. ##

    # Compute the distance to the next recombination event (y).
    y = draw_y(rho=rho, Lx=L0, ploidy=ploidy)
    # While we are still within the sequence intervals.
    while (x + y) < 1:
        # Intialize the new right position in the edge table.
        edge_table.right = np.where(edge_table.right == 1.0, (x + y), edge_table.right)
        # Intialize the edge information for the current tree.
        c_tree_edge_table, c_tree_dicc, Lx = extract_tree_info(
            node_table=node_table, edge_table=edge_table, right=(x + y),
        )

    ## (3) Determine the location (ie what edge), and the age of the recombination event (g). ##

        # Determine g and its location on the current tree.
        rec_edge_key, g = draw_g(c_tree_dicc=c_tree_dicc, Lx=Lx)

    ## (4) Overlay the recombination event at time g and allow the branch below g to coalesce elsewhere on the tree. ##

        # Deteremine the location and time of the recombining coalescent event.
        coal_edge_key, coal_time, root_node = draw_coal_smc(
            node_table=node_table, c_tree_edge_table=c_tree_edge_table, c_tree_dicc=c_tree_dicc,
            rec_edge_key=rec_edge_key, g=g, Ne=Ne, ploidy=ploidy,
        )

    ## (5) Prune the old branch above g and graft the new branch to construct the next tree at position x+y. ##

        # Perform the SPR that is compatible with tskit.
        node_table, edge_table, node_id, _ = ts_spr(
            k=k, node_table=node_table, edge_table=edge_table, c_tree_edge_table=c_tree_edge_table,
            c_tree_dicc=c_tree_dicc, rec_edge_key=rec_edge_key, coal_edge_key=coal_edge_key,
            coal_time=coal_time, right=(x + y), root_node=root_node, node_id=node_id,
        )

    ## (6) Reset the new interval x=x+y, intialize the new tree as the current tree T(x), and compute the compute the total branch length L(x). ##

        # Reset the new left interval (x).
        x = (x + y)
        # Extract the total branch length of the next tree.
        _, _, Lx = extract_tree_info(node_table=node_table, edge_table=edge_table, right=1.0)
        # Compute the distance to the next recombination event (y).
        y = draw_y(rho=rho, Lx=Lx, ploidy=ploidy)
    
    # Validate the new tables to be converted to a tree-sequence.
    smc_tables.sort()
    # Build the new index for the tree-sequence.
    smc_tables.build_index()
    # Create the new tree-sequence.
    smc_ts = smc_tables.tree_sequence()
    # If exporting the tree sequence was specified.
    if export:
        # Export the tree-sequence.
        smc_ts.dump(path)
    return smc_ts

# Define a function to simulate a tree-sequence under the SMC' model.
def sim_ts_smc_prime(k, Ne, rho, ploidy, seed=None, export=False, path='./smc_prime_ts'):
    """
    Returns a tree-sequence simulated using the Sequentially Markovian Coalescent (SMC').
    
    k      -- Number of chromosomes to simulate.
    Ne     -- Effective population size.
    rho    -- Population recombination rate.
    ploidy -- Haploid or diploid coalescent units.
    seed   -- Random seed for simulating T_{0}.
    export -- Do you want to export the tree-sequence?
    path   -- Path to export the tree-sequence. 
    """
    ## (1) Intialize the first tree, T(x)=T_{0}, at position x=0, and compute the total branch length L(x)=L_{0}. ##

    # Intialize the start position.
    x = 0
    # Simulate a tree (T_{0}) under the standard coalescent at point x=0.
    smc_tables, root_id, L0 = sim_T0(k=k, Ne=Ne, ploidy=ploidy, seed=seed)
    # Extract the tables that we will need to edit.
    node_table = smc_tables.nodes
    edge_table = smc_tables.edges
    # Intialize the next node id for tree-seq compatability.
    node_id = root_id + 1
    # Intialize a variable to count the number of hidden recombination events.
    n_hidden = 0
    
    ## (2) Generate the distance, y=exp[(rho/2)L(x)], to the next recombination event. ##

    # Compute the distance to the next recombination event (y).
    y = draw_y(rho=rho, Lx=L0, ploidy=ploidy)
    # While we are still within the sequence intervals.
    while (x + y) < 1:
        # Intialize the new right position in the edge table.
        edge_table.right = np.where(edge_table.right == 1.0, (x + y), edge_table.right)
        # Intialize the edge information for the current tree.
        c_tree_edge_table, c_tree_dicc, Lx = extract_tree_info(
            node_table=node_table, edge_table=edge_table, right=(x + y),
        )

    ## (3) Determine the location (ie what edge), and the age of the recombination event (g). ##

        # Determine g and its location on the current tree.
        rec_edge_key, g = draw_g(c_tree_dicc=c_tree_dicc, Lx=Lx)

    ## (4) Overlay the recombination event at time g and allow the branch below g to coalesce elsewhere on the tree. ##

        # Deteremine the location and time of the recombining coalescent event.
        coal_edge_key, coal_time, root_node = draw_coal_smc_prime(
            node_table=node_table, c_tree_edge_table=c_tree_edge_table, c_tree_dicc=c_tree_dicc,
            g=g, Ne=Ne, ploidy=ploidy,
        )

    ## (5) Prune the old branch above g and graft the new branch to construct the next tree at position x+y. ##

        # Perform the SPR that is compatible with tskit.
        node_table, edge_table, node_id, hidden = ts_spr(
            k=k, node_table=node_table, edge_table=edge_table, c_tree_edge_table=c_tree_edge_table,
            c_tree_dicc=c_tree_dicc, rec_edge_key=rec_edge_key, coal_edge_key=coal_edge_key,
            coal_time=coal_time, right=(x + y), root_node=root_node, node_id=node_id,
        )
        # Update the hidden recombination event variable.
        n_hidden += hidden

    ## (6) Reset the new interval x=x+y, intialize the new tree as the current tree T(x), and compute the compute the total branch length L(x). ##

        # Reset the new left interval (x).
        x = (x + y)
        # Extract the total branch length of the next tree.
        _, _, Lx = extract_tree_info(node_table=node_table, edge_table=edge_table, right=1.0)
        # Compute the distance to the next recombination event (y).
        y = draw_y(rho=rho, Lx=Lx, ploidy=ploidy)
    
    # Validate the new tables to be converted to a tree-sequence.
    smc_tables.sort()
    # Build the new index for the tree-sequence.
    smc_tables.build_index()
    # Create the new tree-sequence.
    smc_ts = smc_tables.tree_sequence()
    # If exporting the tree sequence was specified.
    if export:
        # Export the tree-sequence.
        smc_ts.dump(path)
    return smc_ts, n_hidden