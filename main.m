import "voronoi.m" : initialPerfectForm, neighbour, voronoiAlgorithm, equivalent, automorphisms;
import "symmetricSpace.m" : minimalVectors, toHermitians;
import "init.m" : orderBasis, matrixRationalToImaginary, theta, rationalMatrixEmbedding, matricesB, Dagger, orderBasisCoordinates;
import "cohomology.m" : retract, retractFacets, chainMap, complex, facesEqual, orientationPreserving, orientability;
import "polytope.m" : facetsAsForms, barycentre;

//Voronoi homology

perfectForms := voronoiAlgorithm();
voronoiFacets := [];
for form in perfectForms do
	_, minVecs := minimalVectors(form);
	Append(~voronoiFacets, toHermitians(minVecs));
end for;

print #perfectForms, "classes of perfect forms";

cells, representatives, equivalenceIndices, equivalenceWitnesses := complex(voronoiFacets);
oriented := orientability(representatives);
cells := Reverse(cells);
representatives := Reverse(representatives);
equivalenceIndices := Reverse(equivalenceIndices);
equivalenceWitnesses := Reverse(equivalenceWitnesses);
oriented := Reverse(oriented);

//for cell in representatives do
//	print cell;
//end for;

chainMaps := [* *];
Append(~chainMaps, chainMap(cells, representatives, oriented, equivalenceIndices, equivalenceWitnesses, 1)); //map 2 to 1
Append(~chainMaps, chainMap(cells, representatives, oriented, equivalenceIndices, equivalenceWitnesses, 2)); //map 3 to 2
print chainMaps;
print oriented;

ranks := [];
count := 0;
for j in [1..#oriented[1]] do
	if oriented[1][j] then
		count +:= 1;
	end if;
end for;
kernelDims := [count];

for i in [1..2] do
	Append(~ranks, Rank(chainMaps[i]));
	
	count := 0;
	for j in [1..#oriented[i+1]] do
		if not oriented[i+1][j] then
			count +:= 1;
		end if;
	end for;
	
	Append(~kernelDims, NumberOfRows(KernelMatrix(chainMaps[i]))-count); //3,7,8; 3,4,6 vs 3,8,9; 3,4,7 god wtf
end for;
Append(~ranks, 0);
print ranks, kernelDims;

for i in [1..3] do
	print "dim H _", i-1, "=", kernelDims[i]-ranks[i];
end for;

//---------------Well-rounded retract cohomology------------------------
/*
perfectForms := voronoiAlgorithm();
print #perfectForms, "classes of perfect forms";
_, minVecs := minimalVectors(perfectForms[1]);

cells, representatives, orientability, equivalenceIndices, equivalenceWitnesses := retract(perfectForms);


for i in [1..#cells] do
        print #representatives[i], "cells of dimension", i-1;
        print orientability[i];
        //print equivalenceIndices[i];
end for;


chainMaps := [* *];
Append(~chainMaps, chainMap(cells, representatives, orientability, equivalenceIndices, equivalenceWitnesses, 1));
Append(~chainMaps, chainMap(cells, representatives, orientability, equivalenceIndices, equivalenceWitnesses, 2));
//print chainMaps;
//print chainMaps[2] * chainMaps[1] eq 0;


ranks := [0];
kernelDims := [];
for i in [1..2] do
	Append(~ranks, Rank(chainMaps[i]));
	
	count := 0;
	for j in [1..#orientability[i]] do
		if not orientability[i][j] then
			print "btw there is some non-orientable stuff going on";
			count +:= 1;
		end if;
	end for;
	
	Append(~kernelDims, NumberOfRows(KernelMatrix(Transpose(chainMaps[i])))-count);
end for;

count := 0;
for j in [1..#orientability[3]] do
	if orientability[3][j] then
		count +:= 1;
	end if;
end for;
Append(~kernelDims, count);
print ranks, kernelDims;

for i in [1..#ranks] do
	print "dim H ^", i-1, "=", kernelDims[i]-ranks[i];
end for;*/
