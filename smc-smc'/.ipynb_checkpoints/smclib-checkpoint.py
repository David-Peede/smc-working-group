# Import packages.
import copy
import msprime
import numpy as np

# Intialize a node class.
class Node:
    
    # Intialize the node.
    def __init__(self, node_id, age, node_type, parent=None, l_child=None, r_child=None):
        """
        Node Types
            - 0: leaf node
            - 1: coalescent event node
            - 2: visibile recombination
            - 3: hidden recombination
        """
        self.node_id = node_id
        self.age = age
        self.node_type = node_type
        self.parent = parent
        self.l_child = l_child
        self.r_child = r_child
        self.parent_dist = None
        self.l_child_dist = None
        self.r_child_dist = None
        
    # Define a deep copy method.
    def __deepcopy__(self, memo):
        """
        Return a deepy copy of an isntance of the Node class.
        """
        # Avoid infinite loops
        if id(self) in memo:
            return memo[id(self)]
        # Create a shallow copy of the current node
        copied_node = copy.copy(self)
        memo[id(self)] = copied_node
        # Deep copy children and parent
        copied_node.parent = copy.deepcopy(self.parent, memo)
        copied_node.l_child = copy.deepcopy(self.l_child, memo)
        copied_node.r_child = copy.deepcopy(self.r_child, memo)
        return copied_node
    
    # Define a method to check if a node is a leaf.
    def is_leaf(self):
        """
        True if the node is a leaf, False otherwise.
        """
        return self.node_type == 0
    
    # Define a method to compute the distance to the children.
    def dist_to_children(self):
        """
        Compute the distance from the current node to its children.
        """
        if self.l_child is not None:
            self.l_child_dist = self.age - self.l_child.age
        if self.r_child is not None:
            self.r_child_dist = self.age - self.r_child.age
    
    # Define a method to compute the distance to the parent node.
    def dist_to_parent(self):
        """
        Compute the distance from the current node to its parent.
        """
        if self.parent is not None:
            self.parent_dist = self.parent.age - self.age
            
    # Define a function to initialize distance to parent and children nodes.
    def init_dists(self):
        """
        Intialize the distances to the parent and children nodes.
        """
        self.dist_to_parent()
        self.dist_to_children()


# Intialize a tree class.
class Tree:
    
    # Intialize the tree.
    def __init__(self, left=0.0, right=1.0):
        self.left = left
        self.right = right
        self.root = None
        self.length = None
        self.next_node_id = None
        self.next_rec_id = -1
        self.nodes = {}
        self.edges = {}
        self.upper_bounds = None
        self.recomb_node = None
        self.coal_node = None
        
    def __deepcopy__(self, memo):
        """
        Return a deepy copy of an isntance of the Tree class.
        """
        # Avoid infinite loops.
        if id(self) in memo:
            return memo[id(self)]
        # Create a shallow copy of the tree.
        copied_tree    = copy.copy(self)
        memo[id(self)] = copied_tree
        # Deep copy nodes and edges.
        copied_tree.nodes = copy.deepcopy(self.nodes, memo)
        copied_tree.edges = copy.deepcopy(self.edges, memo)
        # Deep copy unary nodes.
        copied_tree.recomb_node = copy.deepcopy(self.recomb_node, memo)
        copied_tree.coal_node   = copy.deepcopy(self.coal_node, memo)
        return copied_tree
        
    # Define a method to add a node to the tree.
    def add_node(self, node):
        """
        Add a new node to the tree.
        """
        self.nodes[node.node_id] = node
        
    # Define a method to remove a node from the tree.
    def rmv_node(self, node):
        """
        Remove a new node to the tree.
        """
        del self.nodes[node.node_id]
    
    # Define a method to intialize node distances.
    def init_branch_lengths(self):
        """
        Intialize all the branch lengths for the current tree.
        """    
        # For every node.
        for node_id in self.nodes:
            # Intialize branch lengths.
            self.nodes[node_id].init_dists()
        
    # Define a method to intialize the edges on a tree.
    def init_edges(self):
        """
        Intialize all the edges on the current tree.
        """
        # Intialize variables.
        i = 0
        Lx = 0
        upper_bounds = []
        # For every node.
        for node in self.nodes:
            # If the node is not a leaf.
            if not self.nodes[node].is_leaf():
                # Record the interval's upper bound.
                upper_bounds.append(self.nodes[node].age)
                # Intialize the edge for parent -> left child.
                self.edges[i] = {}
                self.edges[i]['parent'] = self.nodes[node].node_id
                self.edges[i]['child'] = self.nodes[node].l_child.node_id
                self.edges[i]['upper'] = self.nodes[node].age
                self.edges[i]['lower'] = self.nodes[node].l_child.age
                self.edges[i]['length'] = self.nodes[node].l_child_dist
                i += 1
                Lx += self.nodes[node].l_child_dist
                # Intialize the edge for parent -> right child.
                self.edges[i] = {}
                self.edges[i]['parent'] = self.nodes[node].node_id
                self.edges[i]['child'] = self.nodes[node].r_child.node_id
                self.edges[i]['upper'] = self.nodes[node].age
                self.edges[i]['lower'] = self.nodes[node].r_child.age
                self.edges[i]['length'] = self.nodes[node].r_child_dist
                i += 1
                Lx += self.nodes[node].r_child_dist
        # Set the tree properties.
        self.upper_bounds = sorted(upper_bounds)
        self.length = Lx
                
    # Define a method to find the root node
    def find_root(self):
        """
        Determine the root node on the current tree.
        """
        root_node = max(self.nodes, key=lambda k: self.nodes[k].age)
        self.root = root_node
        
    # Define a method to reconcile the recoal and associated nodes.
    def perform_spr(self, broken_node, inherited_node, lonely_node, below_node, recoal_node, root_node):
        """
        Perform an SPR operation on the current tree.
        """
        # Intialize the new root.
        new_root = None
    
        ## [1] The broken and below nodes are the root node. ##
        if root_node.node_id == broken_node.node_id and root_node.node_id == below_node.node_id:
            
            # Reconcile the recoal node. #
            # Point the recoal node to its children: {inherited node, lonely node}.
            recoal_node.l_child = inherited_node
            recoal_node.r_child = lonely_node
            # Point the children: {inherited node, lonely node} to the recoal node.
            inherited_node.parent = recoal_node
            lonely_node.parent    = recoal_node
            # Set recoal node as the new root.
            new_root = recoal_node
        
        
        ## [2] The broken node is the below node and the root node is unique. ##
        elif broken_node.node_id == below_node.node_id and broken_node.node_id != root_node.node_id:

            # Reconcile the recoal node's parent node. #
            # Set the parent of the coal node as the parent of the recoal node.
            recoal_node.parent = self.coal_node.parent
            # If the left child of the parent is the broken/below node.
            if recoal_node.parent.l_child.node_id == below_node.node_id:
                # Replace the below node with the recoal node in the children set.
                recoal_node.parent.l_child = recoal_node
            # Else the right child of the parent is the below node.
            else:
                # Replace the below node with the recoal node in the children set.
                recoal_node.parent.r_child = recoal_node
            
            # Reconcile the recoal node. #
            # Point the recoal node to its children: {inherited node, lonely node}.
            recoal_node.l_child = inherited_node
            recoal_node.r_child = lonely_node
            # Point the children: {inherited node, lonely node} to the recoal node.
            inherited_node.parent = recoal_node
            lonely_node.parent    = recoal_node
        
        
        ## [3] The broken node is the root node and the below node is the lonely node. ##
        elif broken_node.node_id == root_node.node_id and below_node.node_id == lonely_node.node_id:

            # Reconcile the recoal node. #
            # Point the recoal node to its children: {inherited node, below node}.
            recoal_node.l_child = inherited_node
            recoal_node.r_child = below_node
            # Point the children: {inherited node, below node} to the recoal node.
            inherited_node.parent = recoal_node
            below_node.parent     = recoal_node
            # Set recoal node as the new root.
            new_root = recoal_node
        
        
        ## [4] The broken node is the root node and the below and lonely nodes are unique. ##
        elif broken_node.node_id == root_node.node_id and below_node.node_id != lonely_node.node_id:

            # Reconcile the recoal node's parent node. #
            # Set the parent of the coal node as the parent of the recoal node.
            recoal_node.parent = self.coal_node.parent
            # If the left child of the parent is the below node.
            if recoal_node.parent.l_child.node_id == below_node.node_id:
                # Replace the below node with the recoal node in the children set.
                recoal_node.parent.l_child = recoal_node
            # Else the right child of the parent is the below node.
            else:
                # Replace the below node with the recoal node in the children set.
                recoal_node.parent.r_child = recoal_node
            
            # Reconcile the recoal node. #
            # Point the recoal node to its children: {inherited node, below node}.
            recoal_node.l_child = inherited_node
            recoal_node.r_child = below_node
            # Point the children: {inherited node, below node} to the recoal node.
            inherited_node.parent = recoal_node
            below_node.parent     = recoal_node
            
            # Reconcile the lonely node. #
            # Set the lonely node as the new root.
            lonely_node.parent = None
            new_root = lonely_node
        
        
        ## [5] The below node is the root node and the broken node is unique. ##
        elif below_node.node_id == root_node.node_id and broken_node.node_id != root_node.node_id:
            
            # Reconcile the lonely node. #
            # Set the parent of the broken node as the new parent of the lonely node.
            lonely_node.parent = broken_node.parent
            # If the left child of the parent is the broken node.
            if lonely_node.parent.l_child.node_id == broken_node.node_id:
                # Replace the broken node with the lonely node in the children set.
                lonely_node.parent.l_child = lonely_node
            # Else the right child of is the broken node.
            else:
                # Replace the broken node with the lonely node in the children set.
                lonely_node.parent.r_child = lonely_node
            
            
            # Reconcile the recoal node. #
            # Point the recoal node to its children: {inherited node, below node}.
            recoal_node.l_child = inherited_node
            recoal_node.r_child = below_node
            # Point the children: {inherited node, below node} to the recoal node.
            inherited_node.parent = recoal_node
            below_node.parent     = recoal_node
            # Set recoal node as the new root.
            new_root = recoal_node
            
          
        ## [6] The broken node, below node, and root node are all unique, and the below node is the lonely node. ##
        elif broken_node.node_id != root_node.node_id and below_node.node_id == lonely_node.node_id:
            
            # Reconcile the recoal node's parent node. #
            # Set the parent of the broken node as the parent of the recoal node.
            recoal_node.parent = broken_node.parent
            # If the left child of the parent is the broken node.
            if recoal_node.parent.l_child.node_id == broken_node.node_id:
                # Replace the broken node with the recoal node in the children set.
                recoal_node.parent.l_child = recoal_node
            # Else the right child of the parent is the broken node.
            else:
                # Replace the broken node with the recoal node in the children set.
                recoal_node.parent.r_child = recoal_node
            
            # Reconcile the recoal node. #
            # Point the recoal node to its children: {inherited node, below node}.
            recoal_node.l_child = inherited_node
            recoal_node.r_child = below_node
            # Point the children: {inherited node, below node} to the recoal node.
            inherited_node.parent = recoal_node
            below_node.parent     = recoal_node
            
        ## [7] All nodes are unique. ##
        else:
            
            # Reconcile the lonely node. #
            # Set the parent of the broken node as the parent of the lonely node.
            lonely_node.parent = broken_node.parent
            # If the left child of is the broken node.
            if lonely_node.parent.l_child.node_id == broken_node.node_id:
                # Replace the broken node with the lonely node in the children set.
                lonely_node.parent.l_child = lonely_node
            # Else the right child of the below node is the broken node.
            else:
                # Replace the broken node with the lonely node in the children set.
                lonely_node.parent.r_child = lonely_node
                
            # Reconcile the recoal node's parent node. #
            # Set the parent of the coal node as the parent of the recoal node.
            recoal_node.parent = self.coal_node.parent
            # If the left child of the parent is the below node.
            if recoal_node.parent.l_child.node_id == below_node.node_id:
                # Replace the below node with the recoal node in the children set.
                recoal_node.parent.l_child = recoal_node
            # Else the right child of the parent is the below node.
            else:
                # Replace the below node with the recoal node in the children set.
                recoal_node.parent.r_child = recoal_node
            
            # Reconcile the recoal node. #
            # Point the recoal node to its children: {inherited node, below node}.
            recoal_node.l_child = inherited_node
            recoal_node.r_child = below_node
            # Point the children: {inherited node, below node} to the recoal node.
            inherited_node.parent = recoal_node
            below_node.parent     = recoal_node
        
        
        ## Tie up the new tree ##
        # Prune the broken node from the tree.
        self.rmv_node(broken_node)
        # Graft the recoal node to the tree.
        self.add_node(recoal_node)
        # Intialize branch lengths for the new tree.
        self.init_branch_lengths()
        # Intialize the edges for the new tree.
        self.init_edges()
        # Intialize the next node ids for the new tree.
        self.init_next_node_id()
        # If we already know the new root node.
        if new_root is not None:
            # Intialize the new root.
            self.root = new_root.node_id
        # Else, we do not know the root.
        else:
            # Intialize the new root.
            self.find_root()

    
    # Define a function to set the next node id.
    def init_next_node_id(self):
        """
        Set the next node id.
        """
        last_coal = self.coal_node
        max_node = max(self.nodes)
        if last_coal is not None:
            self.next_node_id = max([last_coal.node_id, max_node]) + 1
        else:
            self.next_node_id = max_node + 1
    
    # Define a method to recursively construct the newick information.
    def _to_newick_recursive(self, node):
        """
        Recursively construct the Newick information for a given node.
        """
        # Return the leaf id if the node is a leaf.
        if node.is_leaf():
            return f'{node.node_id}'
        # For internal nodes, get the newick information for each child.
        l_child_info = self._to_newick_recursive(node.l_child)
        r_child_info = self._to_newick_recursive(node.r_child)
        return '({}:{},{}:{})'.format(l_child_info, node.l_child_dist, 
                                      r_child_info, node.r_child_dist)
    
    # Define a method to export a tree in newick format.
    def to_newick(self):
        """
        Convert the tree to its Newick format.
        """
        # Start the conversion from the root.
        nwk = self._to_newick_recursive(self.nodes[self.root])
        # The Newick format ends with a semicolon.
        return nwk + ';'
    
# Define a function to intialize a tree from a msprime simulaion.
def init_msp_tree(k, Ne, ploidy, seed=None):
    """
    Returns a Tree object from a msprime simulation.
    
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
    # Intialize the current tree.
    tree = Tree()
    # For ever node.
    for node_id, age in enumerate(ts.tables.nodes.time):
        # If the node is a leaf.
        if age == 0:
            # Intialize the node.
            node = Node(
                node_id=node_id, age=(age / (Ne * ploidy)), node_type=0,
                parent=None, l_child=None, r_child=None,
            )
            # Add the node to the tree.
            tree.add_node(node)
        # Else, the node is an ancestral node.
        else:
            # Intialize the node.
            node = Node(
                node_id=node_id, age=(age / (Ne * ploidy)), node_type=1,
                parent=None, l_child=None, r_child=None,
            )
            # Add the node to the tree.
            tree.add_node(node)
    # For every parent node.
    for parent in np.unique(ts.tables.edges.parent):
        # Find the children of the parent node.
        left_child, right_child = ts.tables.edges[ts.tables.edges.parent == parent].child
        # Update the parent node for the two children.
        tree.nodes[left_child].parent = tree.nodes[parent]
        tree.nodes[right_child].parent = tree.nodes[parent]
        # Update the children nodes for the parent.
        tree.nodes[parent].l_child = tree.nodes[left_child]
        tree.nodes[parent].r_child = tree.nodes[right_child]
    # Intialize branch lengths.
    tree.init_branch_lengths()
    # Intialize the edges for the current tree.
    tree.init_edges()
    # Intialize the root node.
    tree.find_root()
    # Intialize the next node id.
    tree.init_next_node_id()
    # Make a copy of the tree-seq tables for editting.
    ts_tables = ts.dump_tables()
    return tree, ts, ts_tables

# Define a function to determine the distance to the next recombination event.
def draw_y(rho, Lx, ploidy):
    """
    Returns the distance to the next recombination event.
    
    rho    -- Population recombination rate.
    Lx     -- Total branch length of T_{x}.
    ploidy -- Haploid or diploid coalescent units.
    """
    # # Convert the total branch length from generations
    # # to coalescent units.
    # Lx = Lx / ploidy
    # Draw y.
    # y = np.random.exponential((1 / ((rho / ploidy) * Lx)))
    y = np.random.exponential((2 / (rho * Lx)))
    return y

# Define a function to determine the the lineage and age of the recombination event.
def draw_g(tree):
    """
    Returns the recombination event information for the current tree.
    
    tree -- An instance of the current tree.
    """
    # Compute the edge weights (ie edge_length/L(x)).
    edge_weights = [(tree.edges[key]['length'] / tree.length) for key in tree.edges.keys()]
    # Determine which edge will have the recombinatin event.
    rec_edge_key = np.random.choice(list(tree.edges.keys()), p=edge_weights)
    # Determine the age of the recombination event.
    g = np.random.uniform(tree.edges[rec_edge_key]['lower'], tree.edges[rec_edge_key]['upper'])
    return rec_edge_key, g

# Define a function to determine the lineage and age of the next coalescent event for the smc model.
def draw_coal_smc(tree, rec_edge_key, g, Ne, ploidy):
    """
    Returns the edge and coalescent information for the next tree.
    
    tree              -- An instance of the current tree
    rec_edge_key      -- Key of the edge with the recombination event in tree.edges.
    g                 -- Age of the recombination event on the current tree.
    Ne                -- Effective population size.
    ploidy            -- Haploid or diploid coalescent units.
    """
    # Intialize the lower bound of the first interval.
    c_lower_bound = 0
    # Intialize the key of the edge where the coalescent event will occur.
    coal_edge_key = None
    # For every possible coalescent interval.
    for i, c_upper_bound in enumerate(tree.upper_bounds):
        # Determine if the recombination event occurs below the upper bound of the current interval.
        if c_upper_bound > g:
            # Determine the avaiable lineages in this interval.
            available_lineages = [
                key for key in tree.edges.keys() if\
                ((tree.edges[key]['upper'] >= c_upper_bound)\
                & (tree.edges[key]['lower'] <= c_lower_bound)\
                & (key != rec_edge_key))
            ]
            # If there are avaiable lineages.
            if len(available_lineages) > 0:
                # Determine the time of the coalescent event.
                coal_time = g + np.random.exponential((1 / len(available_lineages)))
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
        coal_time = tree.nodes[tree.root].age + np.random.exponential(1)
    return coal_edge_key, coal_time

# Define a function to determine the lineage and age of the next coalescent event for the smc' model.
def draw_coal_smc_prime(tree, g, Ne, ploidy):
    """
    Returns the edge and coalescent information for the next tree.
    
    tree   -- An instance of the current tree
    g      -- Age of the recombination event on the current tree.
    Ne     -- Effective population size.
    ploidy -- Haploid or diploid coalescent units.
    """
    # Intialize the lower bound of the first interval.
    c_lower_bound = 0
    # Intialize the key of the edge where the coalescent event will occur.
    coal_edge_key = None
    # For every possible coalescent interval.
    for i, c_upper_bound in enumerate(tree.upper_bounds):
        # Determine if the recombination event occurs below the upper bound of the current interval.
        if c_upper_bound > g:
            # Determine the avaiable lineages in this interval.
            available_lineages = [
                key for key in tree.edges.keys() if\
                ((tree.edges[key]['upper'] >= c_upper_bound) & (tree.edges[key]['lower'] <= c_lower_bound))
            ] ### YOU CAN ADD THE CONDITION (key != rec_edge_key) FOR SMC ###
            # Determine the time of the coalescent event.
            coal_time = g + np.random.exponential((1 / len(available_lineages)))
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
        coal_time = tree.nodes[tree.root].age + np.random.exponential(1)
    return coal_edge_key, coal_time

# Define a function to simulate an ARG using the SMC.
def sim_smc(k, Ne, rho, ploidy, seed=None):
    """
    Returns a simulated ARG (formatted a tree-sequence dictionary)
    using the Sequentially Markovian Coalescent (SMC).
    
    k      -- Number of chromosomes to simulate.
    Ne     -- Effective population size.
    rho    -- Population recombination rate.
    ploidy -- Haploid or diploid coalescent units.
    seed   -- Random seed for simulating T_{0}.
    export -- Do you want to export the tree-sequence?
    path   -- Path to export the tree-sequence. 
    """
    ## (1) Intialize the first tree, T(x)=T_{0}, at position x=0, and compute the total branch length L(x)=L_{0}. ##

    # Intialize a tree-sequence dictionary.
    ts_dicc = {}
    # Intialize the first tree index.
    tree_idx = 0
    # Intialize the start position.
    x = 0
    # Simulate a tree (T_{0}) under the standard coalescent at point x=0.
    c_tree, ts, ts_tables = init_msp_tree(k=k, Ne=Ne, ploidy=ploidy, seed=seed)

    ## (2) Generate the distance, y=exp[(rho/2)L(x)], to the next recombination event. ##

    # Compute the distance to the next recombination event (y).
    y = draw_y(rho=rho, Lx=c_tree.length, ploidy=ploidy)
    # While we are still within the sequence intervals.
    while (x + y) < 1:
        # Intialize the new right position
        c_tree.right = (x + y)

    ## (3) Determine the location (ie what edge), and the age of the recombination event (g). ##

        # Determine g and its location on the current tree.
        rec_edge_key, g = draw_g(tree=c_tree)

    ## (4) Overlay the recombination event at time g and allow the branch below g to coalesce elsewhere on the tree. ##

        # Deteremine the location and time of the recombining coalescent event.
        coal_edge_key, coal_time = draw_coal_smc(tree=c_tree, rec_edge_key=rec_edge_key, g=g, Ne=Ne, ploidy=ploidy)

    ## (5) Prune the old branch above g and graft the new branch to construct the next tree at position x+y. ##

        ### HIDDEN RECOMBINATION SCENARIO ###
        
        # If the coalescent event is hidden (ie recombination and coalesence occur on the same branch).
        if rec_edge_key == coal_edge_key:

            # Intialize a recombination event node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            recomb_node = Node(
                node_id=c_tree.next_rec_id, age=g, node_type=3,
                parent=c_tree.nodes[c_tree.edges[rec_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[rec_edge_key]['child']], r_child=None,
            )
            # Intialize the coal node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            coal_node = Node(
                node_id=c_tree.next_node_id, age=coal_time, node_type=1,
                parent=c_tree.nodes[c_tree.edges[coal_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[coal_edge_key]['child']], r_child=None,
            )
            # Move the the recombination event node counter back.
            c_tree.next_rec_id -= 1
            # Move the coalescent event node counter forward.
            c_tree.next_node_id += 1
            # Record the recombination event.
            c_tree.recomb_node = recomb_node
            # Record the re-coalesence event.
            c_tree.coal_node = coal_node
            # Add the current tree to the tree-sequence dictionary.
            ts_dicc[tree_idx] = c_tree
            # Move the tree index forward.
            tree_idx += 1
            # Intialize the next tree by copying the current tree.
            n_tree       = copy.deepcopy(c_tree)
            n_tree.left  = (x + y)
            n_tree.right = 1.0
            
        # Else, the coalescent event is not hidden.
        else:

            ### Perform the SPR Algorithim ###

            # Intialize unary nodes—ie the recomb and coal nodes. #
            # Intialize a recombination event node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            recomb_node = Node(
                node_id=c_tree.next_rec_id, age=g, node_type=2,
                parent=c_tree.nodes[c_tree.edges[rec_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[rec_edge_key]['child']], r_child=None,
            )
            # Move the the recombination event node counter back.
            c_tree.next_rec_id -= 1
            # Record the recombination event.
            c_tree.recomb_node = recomb_node
            # If the re-coalescence event is above the root.
            if coal_edge_key == None:
                # Intialize the next node id.
                next_id = c_tree.next_node_id
                # Intialize the coal node for the current tree.
                # Note this is purely for bookkeeping purposes for future conversion to a gARG.
                coal_node = Node(
                    node_id=next_id, age=coal_time, node_type=1,
                    parent=None, l_child=c_tree.nodes[c_tree.root], r_child=None,
                )
                # Record the coal node on the current tree.
                c_tree.coal_node = coal_node
                # Add the current tree to the tree-sequence dictionary.
                ts_dicc[tree_idx] = c_tree
                # Move the tree index forward.
                tree_idx += 1
                # Intialize the next tree by copying the current tree.
                n_tree       = copy.deepcopy(c_tree)
                n_tree.left  = (x + y)
                n_tree.right = 1.0
            # Else the re-coalescence event is on an existing edge.
            else:
                # Intialize the next node id.
                next_id = c_tree.next_node_id
                # Intialize the coal node for the current tree.
                # Note this is purely for bookkeeping purposes for future conversion to a gARG.
                coal_node = Node(
                    node_id=next_id, age=coal_time, node_type=1,
                    parent=c_tree.nodes[c_tree.edges[coal_edge_key]['parent']],
                    l_child=c_tree.nodes[c_tree.edges[coal_edge_key]['child']], r_child=None,
                )
                # Record the coal node on the current tree.
                c_tree.coal_node = coal_node
                # Add the current tree to the tree-sequence dictionary.
                ts_dicc[tree_idx] = c_tree
                # Move the tree index forward.
                tree_idx += 1
                # Intialize the next tree by copying the current tree.
                n_tree       = copy.deepcopy(c_tree)
                n_tree.left  = (x + y)
                n_tree.right = 1.0

            # Intialize nodes of interest. #
            broken_node    = n_tree.recomb_node.parent
            inherited_node = n_tree.recomb_node.l_child
            lonely_node    = broken_node.r_child if broken_node.l_child.node_id == inherited_node.node_id else broken_node.l_child
            below_node     = n_tree.coal_node.l_child
            root_node      = n_tree.nodes[n_tree.root]
            # Intialize the re-coalesence event node to create the next tree.
            recoal_node = Node(
                node_id=next_id, age=coal_time, node_type=1,
                parent=None, l_child=None, r_child=None,
            )
            ## Reconcile the new tree! ##
            n_tree.perform_spr(
                broken_node=broken_node, inherited_node=inherited_node, lonely_node=lonely_node,
                below_node=below_node, recoal_node=recoal_node, root_node=root_node,
            )


        # Set the new tree as the current tree.
        c_tree = n_tree

    ## (6) Reset the new interval x=x+y, intialize the new tree as the current tree T(x), and compute the compute the total branch length L(x). ##

        # Reset the new left interval (x).
        x = (x + y)
        # Compute the distance to the next recombination event (y).
        y = draw_y(rho=rho, Lx=c_tree.length, ploidy=ploidy)
    
    ## Formatting the final tree. ###
    
    # Remove the recombination an re-coalescence nodes from the last tree,
    # that did not expirence recombination.
    c_tree.recomb_node = None
    c_tree.recoal_node = None
    # Add the last tree to the tree-sequence.
    ts_dicc[tree_idx] = c_tree
    return ts_dicc

# Define a function to simulate an ARG using the SMC'.
def sim_smc_prime(k, Ne, rho, ploidy, seed=None):
    """
    Returns a simulated ARG (formatted a tree-sequence dictionary)
    using the Sequentially Markovian Coalescent (SMC').
    
    k      -- Number of chromosomes to simulate.
    Ne     -- Effective population size.
    rho    -- Population recombination rate.
    ploidy -- Haploid or diploid coalescent units.
    seed   -- Random seed for simulating T_{0}.
    export -- Do you want to export the tree-sequence?
    path   -- Path to export the tree-sequence. 
    """
    ## (1) Intialize the first tree, T(x)=T_{0}, at position x=0, and compute the total branch length L(x)=L_{0}. ##

    # Intialize a tree-sequence dictionary.
    ts_dicc = {}
    # Intialize the first tree index.
    tree_idx = 0
    # Intialize the start position.
    x = 0
    # Simulate a tree (T_{0}) under the standard coalescent at point x=0.
    c_tree, ts, ts_tables = init_msp_tree(k=k, Ne=Ne, ploidy=ploidy, seed=seed)

    ## (2) Generate the distance, y=exp[(rho/2)L(x)], to the next recombination event. ##

    # Compute the distance to the next recombination event (y).
    y = draw_y(rho=rho, Lx=c_tree.length, ploidy=ploidy)
    # While we are still within the sequence intervals.
    while (x + y) < 1:
        # Intialize the new right position
        c_tree.right = (x + y)

    ## (3) Determine the location (ie what edge), and the age of the recombination event (g). ##

        # Determine g and its location on the current tree.
        rec_edge_key, g = draw_g(tree=c_tree)

    ## (4) Overlay the recombination event at time g and allow the branch below g to coalesce elsewhere on the tree. ##

        # Deteremine the location and time of the recombining coalescent event.
        coal_edge_key, coal_time = draw_coal_smc_prime(tree=c_tree, g=g, Ne=Ne, ploidy=ploidy)

    ## (5) Prune the old branch above g and graft the new branch to construct the next tree at position x+y. ##

        ### HIDDEN RECOMBINATION SCENARIO ###
        
        # If the coalescent event is hidden (ie recombination and coalesence occur on the same branch).
        if rec_edge_key == coal_edge_key:

            # Intialize a recombination event node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            recomb_node = Node(
                node_id=c_tree.next_rec_id, age=g, node_type=3,
                parent=c_tree.nodes[c_tree.edges[rec_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[rec_edge_key]['child']], r_child=None,
            )
            # Intialize the coal node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            coal_node = Node(
                node_id=c_tree.next_node_id, age=coal_time, node_type=1,
                parent=c_tree.nodes[c_tree.edges[coal_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[coal_edge_key]['child']], r_child=None,
            )
            # Move the the recombination event node counter back.
            c_tree.next_rec_id -= 1
            # Move the coalescent event node counter forward.
            c_tree.next_node_id += 1
            # Record the recombination event.
            c_tree.recomb_node = recomb_node
            # Record the re-coalesence event.
            c_tree.coal_node = coal_node
            # Add the current tree to the tree-sequence dictionary.
            ts_dicc[tree_idx] = c_tree
            # Move the tree index forward.
            tree_idx += 1
            # Intialize the next tree by copying the current tree.
            n_tree       = copy.deepcopy(c_tree)
            n_tree.left  = (x + y)
            n_tree.right = 1.0
            
        # Else, the coalescent event is not hidden.
        else:

            ### Perform the SPR Algorithim ###

            # Intialize unary nodes—ie the recomb and coal nodes. #
            # Intialize a recombination event node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            recomb_node = Node(
                node_id=c_tree.next_rec_id, age=g, node_type=2,
                parent=c_tree.nodes[c_tree.edges[rec_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[rec_edge_key]['child']], r_child=None,
            )
            # Move the the recombination event node counter back.
            c_tree.next_rec_id -= 1
            # Record the recombination event.
            c_tree.recomb_node = recomb_node
            # If the re-coalescence event is above the root.
            if coal_edge_key == None:
                # Intialize the next node id.
                next_id = c_tree.next_node_id
                # Intialize the coal node for the current tree.
                # Note this is purely for bookkeeping purposes for future conversion to a gARG.
                coal_node = Node(
                    node_id=next_id, age=coal_time, node_type=1,
                    parent=None, l_child=c_tree.nodes[c_tree.root], r_child=None,
                )
                # Record the coal node on the current tree.
                c_tree.coal_node = coal_node
                # Add the current tree to the tree-sequence dictionary.
                ts_dicc[tree_idx] = c_tree
                # Move the tree index forward.
                tree_idx += 1
                # Intialize the next tree by copying the current tree.
                n_tree       = copy.deepcopy(c_tree)
                n_tree.left  = (x + y)
                n_tree.right = 1.0
            # Else the re-coalescence event is on an existing edge.
            else:
                # Intialize the next node id.
                next_id = c_tree.next_node_id
                # Intialize the coal node for the current tree.
                # Note this is purely for bookkeeping purposes for future conversion to a gARG.
                coal_node = Node(
                    node_id=next_id, age=coal_time, node_type=1,
                    parent=c_tree.nodes[c_tree.edges[coal_edge_key]['parent']],
                    l_child=c_tree.nodes[c_tree.edges[coal_edge_key]['child']], r_child=None,
                )
                # Record the coal node on the current tree.
                c_tree.coal_node = coal_node
                # Add the current tree to the tree-sequence dictionary.
                ts_dicc[tree_idx] = c_tree
                # Move the tree index forward.
                tree_idx += 1
                # Intialize the next tree by copying the current tree.
                n_tree       = copy.deepcopy(c_tree)
                n_tree.left  = (x + y)
                n_tree.right = 1.0

            # Intialize nodes of interest. #
            broken_node    = n_tree.recomb_node.parent
            inherited_node = n_tree.recomb_node.l_child
            lonely_node    = broken_node.r_child if broken_node.l_child.node_id == inherited_node.node_id else broken_node.l_child
            below_node     = n_tree.coal_node.l_child
            root_node      = n_tree.nodes[n_tree.root]
            # Intialize the re-coalesence event node to create the next tree.
            recoal_node = Node(
                node_id=next_id, age=coal_time, node_type=1,
                parent=None, l_child=None, r_child=None,
            )
            ## Reconcile the new tree! ##
            n_tree.perform_spr(
                broken_node=broken_node, inherited_node=inherited_node, lonely_node=lonely_node,
                below_node=below_node, recoal_node=recoal_node, root_node=root_node,
            )


        # Set the new tree as the current tree.
        c_tree = n_tree

    ## (6) Reset the new interval x=x+y, intialize the new tree as the current tree T(x), and compute the compute the total branch length L(x). ##

        # Reset the new left interval (x).
        x = (x + y)
        # Compute the distance to the next recombination event (y).
        y = draw_y(rho=rho, Lx=c_tree.length, ploidy=ploidy)
    
    ## Formatting the final tree. ###
    
    # Remove the recombination an re-coalescence nodes from the last tree,
    # that did not expirence recombination.
    c_tree.recomb_node = None
    c_tree.recoal_node = None
    # Add the last tree to the tree-sequence.
    ts_dicc[tree_idx] = c_tree
    return ts_dicc

# Define a function to simulate an ARG using the SMC.
def ith_tree_smc(k, Ne, rho, ploidy, seed=None):
    """
    Returns the first 10 trees and the last tree simulated
    using the Sequentially Markovian Coalescent (SMC).
    
    k      -- Number of chromosomes to simulate.
    Ne     -- Effective population size.
    rho    -- Population recombination rate.
    ploidy -- Haploid or diploid coalescent units.
    seed   -- Random seed for simulating T_{0}.
    export -- Do you want to export the tree-sequence?
    path   -- Path to export the tree-sequence. 
    """
    ## (1) Intialize the first tree, T(x)=T_{0}, at position x=0, and compute the total branch length L(x)=L_{0}. ##

    # Intialize a tree-sequence dictionary.
    ts_dicc = {}
    # Intialize the first tree index.
    tree_idx = 0
    # Intialize the start position.
    x = 0
    # Simulate a tree (T_{0}) under the standard coalescent at point x=0.
    c_tree, ts, ts_tables = init_msp_tree(k=k, Ne=Ne, ploidy=ploidy, seed=seed)

    ## (2) Generate the distance, y=exp[(rho/2)L(x)], to the next recombination event. ##

    # Compute the distance to the next recombination event (y).
    y = draw_y(rho=rho, Lx=c_tree.length, ploidy=ploidy)
    # While we are still within the sequence intervals.
    while (x + y) < 1:
        # Intialize the new right position
        c_tree.right = (x + y)

    ## (3) Determine the location (ie what edge), and the age of the recombination event (g). ##

        # Determine g and its location on the current tree.
        rec_edge_key, g = draw_g(tree=c_tree)

    ## (4) Overlay the recombination event at time g and allow the branch below g to coalesce elsewhere on the tree. ##

        # Deteremine the location and time of the recombining coalescent event.
        coal_edge_key, coal_time = draw_coal_smc(tree=c_tree, rec_edge_key=rec_edge_key, g=g, Ne=Ne, ploidy=ploidy)

    ## (5) Prune the old branch above g and graft the new branch to construct the next tree at position x+y. ##

        ### HIDDEN RECOMBINATION SCENARIO ###
        
        # If the coalescent event is hidden (ie recombination and coalesence occur on the same branch).
        if rec_edge_key == coal_edge_key:

            # Intialize a recombination event node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            recomb_node = Node(
                node_id=c_tree.next_rec_id, age=g, node_type=3,
                parent=c_tree.nodes[c_tree.edges[rec_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[rec_edge_key]['child']], r_child=None,
            )
            # Intialize the coal node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            coal_node = Node(
                node_id=c_tree.next_node_id, age=coal_time, node_type=1,
                parent=c_tree.nodes[c_tree.edges[coal_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[coal_edge_key]['child']], r_child=None,
            )
            # Move the the recombination event node counter back.
            c_tree.next_rec_id -= 1
            # Move the coalescent event node counter forward.
            c_tree.next_node_id += 1
            # Record the recombination event.
            c_tree.recomb_node = recomb_node
            # Record the re-coalesence event.
            c_tree.coal_node = coal_node
            # If the current tree is the 11th tree or greater.
            if tree_idx >= 10:
                # Move the tree index forward.
                tree_idx += 1
            # Else, this is one of the first 10 trees.
            else:
                # Add the current tree to the tree-sequence dictionary.
                ts_dicc[tree_idx] = c_tree
                # Move the tree index forward.
                tree_idx += 1
            # Intialize the next tree by copying the current tree.
            n_tree       = copy.deepcopy(c_tree)
            n_tree.left  = (x + y)
            n_tree.right = 1.0
            
        # Else, the coalescent event is not hidden.
        else:

            ### Perform the SPR Algorithim ###

            # Intialize unary nodes—ie the recomb and coal nodes. #
            # Intialize a recombination event node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            recomb_node = Node(
                node_id=c_tree.next_rec_id, age=g, node_type=2,
                parent=c_tree.nodes[c_tree.edges[rec_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[rec_edge_key]['child']], r_child=None,
            )
            # Move the the recombination event node counter back.
            c_tree.next_rec_id -= 1
            # Record the recombination event.
            c_tree.recomb_node = recomb_node
            # If the re-coalescence event is above the root.
            if coal_edge_key == None:
                # Intialize the next node id.
                next_id = c_tree.next_node_id
                # Intialize the coal node for the current tree.
                # Note this is purely for bookkeeping purposes for future conversion to a gARG.
                coal_node = Node(
                    node_id=next_id, age=coal_time, node_type=1,
                    parent=None, l_child=c_tree.nodes[c_tree.root], r_child=None,
                )
                # Record the coal node on the current tree.
                c_tree.coal_node = coal_node
                # If the current tree is the 11th tree or greater.
                if tree_idx >= 10:
                    # Move the tree index forward.
                    tree_idx += 1
                # Else, this is one of the first 10 trees.
                else:
                    # Add the current tree to the tree-sequence dictionary.
                    ts_dicc[tree_idx] = c_tree
                    # Move the tree index forward.
                    tree_idx += 1
                # Intialize the next tree by copying the current tree.
                n_tree       = copy.deepcopy(c_tree)
                n_tree.left  = (x + y)
                n_tree.right = 1.0
            # Else the re-coalescence event is on an existing edge.
            else:
                # Intialize the next node id.
                next_id = c_tree.next_node_id
                # Intialize the coal node for the current tree.
                # Note this is purely for bookkeeping purposes for future conversion to a gARG.
                coal_node = Node(
                    node_id=next_id, age=coal_time, node_type=1,
                    parent=c_tree.nodes[c_tree.edges[coal_edge_key]['parent']],
                    l_child=c_tree.nodes[c_tree.edges[coal_edge_key]['child']], r_child=None,
                )
                # Record the coal node on the current tree.
                c_tree.coal_node = coal_node
                # If the current tree is the 11th tree or greater.
                if tree_idx >= 10:
                    # Move the tree index forward.
                    tree_idx += 1
                # Else, this is one of the first 10 trees.
                else:
                    # Add the current tree to the tree-sequence dictionary.
                    ts_dicc[tree_idx] = c_tree
                    # Move the tree index forward.
                    tree_idx += 1
                # Intialize the next tree by copying the current tree.
                n_tree       = copy.deepcopy(c_tree)
                n_tree.left  = (x + y)
                n_tree.right = 1.0

            # Intialize nodes of interest. #
            broken_node    = n_tree.recomb_node.parent
            inherited_node = n_tree.recomb_node.l_child
            lonely_node    = broken_node.r_child if broken_node.l_child.node_id == inherited_node.node_id else broken_node.l_child
            below_node     = n_tree.coal_node.l_child
            root_node      = n_tree.nodes[n_tree.root]
            # Intialize the re-coalesence event node to create the next tree.
            recoal_node = Node(
                node_id=next_id, age=coal_time, node_type=1,
                parent=None, l_child=None, r_child=None,
            )
            ## Reconcile the new tree! ##
            n_tree.perform_spr(
                broken_node=broken_node, inherited_node=inherited_node, lonely_node=lonely_node,
                below_node=below_node, recoal_node=recoal_node, root_node=root_node,
            )


        # Set the new tree as the current tree.
        c_tree = n_tree

    ## (6) Reset the new interval x=x+y, intialize the new tree as the current tree T(x), and compute the compute the total branch length L(x). ##

        # Reset the new left interval (x).
        x = (x + y)
        # Compute the distance to the next recombination event (y).
        y = draw_y(rho=rho, Lx=c_tree.length, ploidy=ploidy)
    
    ## Formatting the final tree. ###
    
    # Remove the recombination an re-coalescence nodes from the last tree,
    # that did not expirence recombination.
    c_tree.recomb_node = None
    c_tree.recoal_node = None
    # Add the last tree to the tree-sequence.
    ts_dicc[tree_idx] = c_tree
    return ts_dicc

# Define a function to simulate an ARG using the SMC'.
def ith_tree_smc_prime(k, Ne, rho, ploidy, seed=None):
    """
    Returns the first 10 trees and the last tree simulated
    using the Sequentially Markovian Coalescent (SMC').
    
    k      -- Number of chromosomes to simulate.
    Ne     -- Effective population size.
    rho    -- Population recombination rate.
    ploidy -- Haploid or diploid coalescent units.
    seed   -- Random seed for simulating T_{0}.
    export -- Do you want to export the tree-sequence?
    path   -- Path to export the tree-sequence. 
    """
    ## (1) Intialize the first tree, T(x)=T_{0}, at position x=0, and compute the total branch length L(x)=L_{0}. ##

    # Intialize a tree-sequence dictionary.
    ts_dicc = {}
    # Intialize the first tree index.
    tree_idx = 0
    # Intialize the start position.
    x = 0
    # Simulate a tree (T_{0}) under the standard coalescent at point x=0.
    c_tree, ts, ts_tables = init_msp_tree(k=k, Ne=Ne, ploidy=ploidy, seed=seed)

    ## (2) Generate the distance, y=exp[(rho/2)L(x)], to the next recombination event. ##

    # Compute the distance to the next recombination event (y).
    y = draw_y(rho=rho, Lx=c_tree.length, ploidy=ploidy)
    # While we are still within the sequence intervals.
    while (x + y) < 1:
        # Intialize the new right position
        c_tree.right = (x + y)

    ## (3) Determine the location (ie what edge), and the age of the recombination event (g). ##

        # Determine g and its location on the current tree.
        rec_edge_key, g = draw_g(tree=c_tree)

    ## (4) Overlay the recombination event at time g and allow the branch below g to coalesce elsewhere on the tree. ##

        # Deteremine the location and time of the recombining coalescent event.
        coal_edge_key, coal_time = draw_coal_smc_prime(tree=c_tree, g=g, Ne=Ne, ploidy=ploidy)

    ## (5) Prune the old branch above g and graft the new branch to construct the next tree at position x+y. ##

        ### HIDDEN RECOMBINATION SCENARIO ###
        
        # If the coalescent event is hidden (ie recombination and coalesence occur on the same branch).
        if rec_edge_key == coal_edge_key:

            # Intialize a recombination event node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            recomb_node = Node(
                node_id=c_tree.next_rec_id, age=g, node_type=3,
                parent=c_tree.nodes[c_tree.edges[rec_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[rec_edge_key]['child']], r_child=None,
            )
            # Intialize the coal node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            coal_node = Node(
                node_id=c_tree.next_node_id, age=coal_time, node_type=1,
                parent=c_tree.nodes[c_tree.edges[coal_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[coal_edge_key]['child']], r_child=None,
            )
            # Move the the recombination event node counter back.
            c_tree.next_rec_id -= 1
            # Move the coalescent event node counter forward.
            c_tree.next_node_id += 1
            # Record the recombination event.
            c_tree.recomb_node = recomb_node
            # Record the re-coalesence event.
            c_tree.coal_node = coal_node
            # If the current tree is the 11th tree or greater.
            if tree_idx >= 10:
                # Move the tree index forward.
                tree_idx += 1
            # Else, this is one of the first 10 trees.
            else:
                # Add the current tree to the tree-sequence dictionary.
                ts_dicc[tree_idx] = c_tree
                # Move the tree index forward.
                tree_idx += 1
            # Intialize the next tree by copying the current tree.
            n_tree       = copy.deepcopy(c_tree)
            n_tree.left  = (x + y)
            n_tree.right = 1.0
            
        # Else, the coalescent event is not hidden.
        else:

            ### Perform the SPR Algorithim ###

            # Intialize unary nodes—ie the recomb and coal nodes. #
            # Intialize a recombination event node for the current tree.
            # Note this is purely for bookkeeping purposes for future conversion to a gARG.
            recomb_node = Node(
                node_id=c_tree.next_rec_id, age=g, node_type=2,
                parent=c_tree.nodes[c_tree.edges[rec_edge_key]['parent']],
                l_child=c_tree.nodes[c_tree.edges[rec_edge_key]['child']], r_child=None,
            )
            # Move the the recombination event node counter back.
            c_tree.next_rec_id -= 1
            # Record the recombination event.
            c_tree.recomb_node = recomb_node
            # If the re-coalescence event is above the root.
            if coal_edge_key == None:
                # Intialize the next node id.
                next_id = c_tree.next_node_id
                # Intialize the coal node for the current tree.
                # Note this is purely for bookkeeping purposes for future conversion to a gARG.
                coal_node = Node(
                    node_id=next_id, age=coal_time, node_type=1,
                    parent=None, l_child=c_tree.nodes[c_tree.root], r_child=None,
                )
                # Record the coal node on the current tree.
                c_tree.coal_node = coal_node
                # If the current tree is the 11th tree or greater.
                if tree_idx >= 10:
                    # Move the tree index forward.
                    tree_idx += 1
                # Else, this is one of the first 10 trees.
                else:
                    # Add the current tree to the tree-sequence dictionary.
                    ts_dicc[tree_idx] = c_tree
                    # Move the tree index forward.
                    tree_idx += 1
                # Intialize the next tree by copying the current tree.
                n_tree       = copy.deepcopy(c_tree)
                n_tree.left  = (x + y)
                n_tree.right = 1.0
            # Else the re-coalescence event is on an existing edge.
            else:
                # Intialize the next node id.
                next_id = c_tree.next_node_id
                # Intialize the coal node for the current tree.
                # Note this is purely for bookkeeping purposes for future conversion to a gARG.
                coal_node = Node(
                    node_id=next_id, age=coal_time, node_type=1,
                    parent=c_tree.nodes[c_tree.edges[coal_edge_key]['parent']],
                    l_child=c_tree.nodes[c_tree.edges[coal_edge_key]['child']], r_child=None,
                )
                # Record the coal node on the current tree.
                c_tree.coal_node = coal_node
                # If the current tree is the 11th tree or greater.
                if tree_idx >= 10:
                    # Move the tree index forward.
                    tree_idx += 1
                # Else, this is one of the first 10 trees.
                else:
                    # Add the current tree to the tree-sequence dictionary.
                    ts_dicc[tree_idx] = c_tree
                    # Move the tree index forward.
                    tree_idx += 1
                # Intialize the next tree by copying the current tree.
                n_tree       = copy.deepcopy(c_tree)
                n_tree.left  = (x + y)
                n_tree.right = 1.0

            # Intialize nodes of interest. #
            broken_node    = n_tree.recomb_node.parent
            inherited_node = n_tree.recomb_node.l_child
            lonely_node    = broken_node.r_child if broken_node.l_child.node_id == inherited_node.node_id else broken_node.l_child
            below_node     = n_tree.coal_node.l_child
            root_node      = n_tree.nodes[n_tree.root]
            # Intialize the re-coalesence event node to create the next tree.
            recoal_node = Node(
                node_id=next_id, age=coal_time, node_type=1,
                parent=None, l_child=None, r_child=None,
            )
            ## Reconcile the new tree! ##
            n_tree.perform_spr(
                broken_node=broken_node, inherited_node=inherited_node, lonely_node=lonely_node,
                below_node=below_node, recoal_node=recoal_node, root_node=root_node,
            )


        # Set the new tree as the current tree.
        c_tree = n_tree

    ## (6) Reset the new interval x=x+y, intialize the new tree as the current tree T(x), and compute the compute the total branch length L(x). ##

        # Reset the new left interval (x).
        x = (x + y)
        # Compute the distance to the next recombination event (y).
        y = draw_y(rho=rho, Lx=c_tree.length, ploidy=ploidy)
    
    ## Formatting the final tree. ###
    
    # Remove the recombination an re-coalescence nodes from the last tree,
    # that did not expirence recombination.
    c_tree.recomb_node = None
    c_tree.recoal_node = None
    # Add the last tree to the tree-sequence.
    ts_dicc[tree_idx] = c_tree
    return ts_dicc

# Define a function to extract the number of recombination events from a tskit tree-sequence gARG.
def R_g_arg(ts):
    return np.unique(ts.tables.nodes[ts.tables.nodes.flags > 1].time).size

# Define a function to extract the number of recombination events from a tree-sequence ARG dictionary.
def R_arg_dicc(ts_dicc):
    return len(ts_dicc) - 1

# Define a function to extract the tmrca from the ith tree in a tskit tree-sequence.
def ith_tmrca_ts(ts, ith_tree):
    # If the ith tree exists.
    if ts.num_trees > ith_tree:
        # Extract the tmrca.
        ith_tmrca = ts.at_index(ith_tree).time(ts.at_index(ith_tree).root)
        # Find the width of the tree.
        ith_width = ts.at_index(ith_tree).interval.right - ts.at_index(ith_tree).interval.left
    # Else, the ith tree doesn't exist.
    else:
        ith_tmrca = np.nan
        ith_width = np.nan
    return ith_tmrca, ith_width

# Define a function to extract the tmrca from the ith tree in a tree-sequence dictionary.
def ith_tmrca_ts_dicc(ts_dicc, ith_tree):
    # If the ith tree is the last tree.
    if ith_tree == -1:
        ith_tmrca = ts_dicc[max(ts_dicc)].nodes[ts_dicc[max(ts_dicc)].root].age
        # Find the width of the tree.
        ith_width = ts_dicc[max(ts_dicc)].right - ts_dicc[max(ts_dicc)].left
    # Else-if the ith tree exists.
    elif len(ts_dicc) > ith_tree:
        # Extract the tmrca.
        ith_tmrca = ts_dicc[ith_tree].nodes[ts_dicc[ith_tree].root].age
        # Find the width of the tree.
        ith_width = ts_dicc[ith_tree].right - ts_dicc[ith_tree].left
    # Else, the ith tree doesn't exist.
    else:
        ith_tmrca = np.nan
        ith_width = np.nan
    return ith_tmrca, ith_width