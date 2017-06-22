============================================
Geoemtry Processing Introduction
Solution by Jonathan Lesser
Submitted January 29,2017
============================================

---------
edges.cpp
---------
To generate the matrix E, I form every possible edge as a 
std::pair<int,int> by looking at each row of the matrix F
and pairing off the indices of the adjacent vertices, which
we know since the vertices determine the orientation of the
face through their order. I make sure that the first element
in the pair is less than the second, to allow for well-defined
comparison of edges. Then I insert the edge pair into a set 
of pairs and as a row in E if the edge pair is not already in 
the set of pairs, to guarantee there will be no duplicate
edge in E.

Now that I have finished my implementation and I see a
decent amount of simmilar looking code, I am thinking that
there might be a prettier (and possibly more efficient) way
to do this task by converting the matrix to an array and
finding a clever way to select the right vertices to form 
the edges, but I have not worked through this.

------------------------
euler_characteristic.cpp
------------------------
The Euler characteristic Chi is calculated from the number
of faces, edges, and vertices by a simple equation, so the
solution follows quickly, since the number of faces is equal
to the number of rows of F, the number of edges is equal to
the number of rows of edges(F), and the number of vertices
is equal to one plus the maximum integer contained in F.
The maximum integer i contained in F is the index of the last
vertex in the matrix of vertices, so the number of vertices
is i+1.

-------------------------------
Difficulties and clarifications
-------------------------------

My primary difficulties stemmed from rustiness with C++,
so I anticipate they will wear off and I will be better
on the next assignment. I do wonder how important time
complexity is and how focused we should be on finding
very fast algorithms? The assignment was very clear in
its instructions.
