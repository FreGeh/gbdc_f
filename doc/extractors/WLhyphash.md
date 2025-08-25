1-WL on the incidence hypergraph:
* edge label = hash(kind, size, multiset of incident vertex labels)
* vertex label = hash(polarity_bit, degree, multiset of incident edge labels)

Sorting the neighbor labels makes the signature a multiset, which is isomorphism-invariant.