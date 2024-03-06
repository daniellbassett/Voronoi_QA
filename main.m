import "voronoi.m" : initialPerfectForm, neighbour, voronoiAlgorithm, equivalent, automorphisms;
import "symmetricSpace.m" : minimalVectors, toHermitians;
import "init.m" : orderBasis, matrixRationalToQuaternion, rationalMatrixEmbedding, matricesB, Dagger, orderBasisCoordinates;
import "cohomology.m" : retract, retractFacets, chainMap;
import "polytope.m" : facetsAsForms, barycentre;


perfectForms := voronoiAlgorithm();
print #perfectForms, "classes of perfect forms";
_, minVecs := minimalVectors(perfectForms[1]);

cells, representatives, orientability, equivalenceIndices, equivalenceWitnesses := retract(perfectForms);

for i in [1..#cells] do
        print #representatives[i], "cells of dimension", i-1;
        print orientability[i];
        print equivalenceIndices[i];
end for;

chainMaps := [* *];

ranks := [0];
kernelDims := [];
for i in [1..4] do
	Append(~chainMaps, chainMap(cells, representatives, orientability, equivalenceIndices, equivalenceWitnesses, i));
	Append(~ranks, Rank(chainMaps[i]));
	
	count := 0;
	for j in [1..#orientability[i]] do
		if not orientability[i][j] then
			count +:= 1;
		end if;
	end for;
	
	Append(~kernelDims, NumberOfRows(KernelMatrix(Transpose(chainMaps[i])))-count);
end for;

count := 0;
for j in [1..#orientability[5]] do
	if orientability[5][j] then
		count +:= 1;
	end if;
end for;
Append(~kernelDims, count);
print ranks, kernelDims;
