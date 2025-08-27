1-WL on the incidence hypergraph:
* edge label = hash(kind, size, multiset of incident vertex labels)
* vertex label = hash(polarity_bit, degree, multiset of incident edge labels)

Sorting the neighbor labels makes the signature a multiset, which is isomorphism-invariant.



possible variations:
what information goes into the hashing:
- counts of distinct colours each round
- sorted histogram of class sizes
- (id, count) pairs per vertex, edge_varpair, edge_clause (for every iteration / only last)
- full sorted list of all ids