import "init.m" : hermDim, matricesB, Dagger, matrixRationalToImaginary;
import "polytope.m" : facetsAsForms, barycentre, toPolytope;
import "symmetricSpace.m" : toHermitians, minimalVectors, formBasis, linearCombinations;
import "voronoi.m" : equivalent, automorphisms;

function facesEqual(face1, face2)
	if #face1 eq #face2 then
		for vertex in face1 do
			if not vertex in face2 then
				return false;
			end if;
		end for;
		
		for vertex in face2 do
			if not vertex in face1 then
				return false;
			end if;
		end for;
	else
		return false;
	end if;
	
	return true;
end function;

function retractFacets(perfectForms) //assumes for now that all minimal vectors are primitive wrt the lattice
	//Facets of well-rounded retract of a cone correspond to flags of its faces, with the final 0-dimensional face removed
	
	//List of faces of cones
	//eventually: faces := [[[cone1facet], [cone1codim1_1, cone1codim1_2], ...], [[cone2facet], ...]]
	faces := [];
	//cones
	for form in perfectForms do
		_, vecs := minimalVectors(form);
		Append(~faces, [[toHermitians(vecs)]]);
	end for;
	
	//Lower dimensional faces
	for i in [1..#perfectForms] do
		for j in [2..hermDim-1] do
			Append(~faces[i], []);
			
			for face in faces[i][j-1] do //All dimension 1 higher faces in the ith cone
				facets := facetsAsForms(face);
				for facet in facets do
					new := true;
					for previousFacet in faces[i][j] do
						if facesEqual(facet, previousFacet) then
							new := false;
							break;
						end if;
					end for;
					
					if new then
						Append(~faces[i][j], facet);
					end if;
				end for;
			end for;
		end for;
	end for;
	
	//Flag construction
	//eventually: flags := [[[[cone1facet]], [[cone1facet,cone1codim1_1],[cone1facet,cone1codim1_2]]], []]
	flags := [[[faces[i][1]]] : i in [1..#perfectForms]];
	for i in [1..#perfectForms] do
		for j in [2..hermDim-1] do //recursively construct flags of length j within the ith cone
			Append(~flags[i], []);
			
			for sequence in flags[i][j-1] do
				facets := facetsAsForms(sequence[j-1]);
				for facet in facets do
					Append(~flags[i][j], Append(sequence, facet));
				end for;
			end for;
			//print i, j, #flags[i][j];
		end for;
	end for;
	
	//Convert flags into cells by association with barycentres
	flags := [flags[i][hermDim-1] : i in [1..#perfectForms]];
	
	facets := [];
	for i in [1..#perfectForms] do
		for j in [1..#flags[i]] do
			Append(~facets, [barycentre(flags[i][j][k]) : k in [1..hermDim-1]]);
		end for;
	end for;
	
	//Find equivalence classes of the above cells
	facetClasses := [];
	count := 1;
	for facet in facets do
		count +:= 1;
		new := true;
		for previousFacet in facetClasses do
			if equivalent(barycentre(facet), barycentre(previousFacet)) then
				new := false;
				break;
			end if;
		end for;
		
		if new then
			Append(~facetClasses, facet);
		end if;
	end for;
	
	return facetClasses;
end function;

function complex(facets) //Generates complex from its facets
	cells := [facets];
	
	representatives := cells;
	equivalenceIndices := [[1..#facets]];
	equivalenceWitnesses := [[matricesB ! 1 : i in [1..#facets]]];
	
	for i in [2..hermDim-1] do
		//Generate codim i cells from codim i-1 cells
		Append(~cells, []);
		
		for cell in representatives[i-1] do
			cells[i] := cells[i] cat facetsAsForms(cell);
		end for;
		
		//Form classes of codim i cells
		Append(~representatives, []);
		Append(~equivalenceIndices, []);
		Append(~equivalenceWitnesses, []);
		
		for cell in cells[i] do
			bary := barycentre(cell);
			new := true;
			
			for k in [1..#representatives[i]] do
				equiv, conjugate := equivalent(barycentre(representatives[i][k]), bary); //conjugate takes bary to barycentre(representatives[i][k])
				
				if equiv then
					conjugate := matrixRationalToImaginary(conjugate);
					
					Append(~equivalenceIndices[i], k);
					Append(~equivalenceWitnesses[i], conjugate);
					
					new := false;
					break;
				end if;
			end for;
			
			if new then
				Append(~representatives[i], cell);
				Append(~equivalenceIndices[i], #representatives[i]);
				Append(~equivalenceWitnesses[i], matricesB ! 1);
			end if;
		end for;
		
		print "Calculated representatives in codimension", i-1;
	end for;
	
	return cells, representatives, equivalenceIndices, equivalenceWitnesses;
end function;

function orientationPreserving(G, S) //Determines if G preserves the orientation on the space spanned by S
	basis := formBasis(S);
	
	for gamma in G do
		action := linearCombinations(basis, [gamma * minForm * Dagger(gamma) : minForm in basis]);
		
		det := Determinant(MatrixRing(Rationals(), #basis) ! action);
		if det lt 0 then
			return false;
		end if;
	end for;
	
	return true;
end function;

function orientability(representatives)
	orientations := [];
	for i in [1..#representatives] do
		Append(~orientations, []);
		
		for representative in representatives[i] do
			G := automorphisms(barycentre(representative));
			orientable := orientationPreserving(G, representative);
			Append(~orientations[i], orientable);
		end for;
	end for;
	
	return orientations;
end function;

function retract(perfectForms)
	facets := retractFacets(perfectForms);
	print "Top dimensional facets of retract calculated";
	cells, representatives, equivalenceIndices, equivalenceWitnesses := complex(facets);
	oriented := orientability(representatives);
	
	return Reverse(cells), Reverse(representatives), Reverse(oriented), Reverse(equivalenceIndices), Reverse(equivalenceWitnesses);
end function;

function chainMap(cells, representatives, orientability, equivalenceIndices, equivalenceWitnesses, lowDim) //Calculates the orientation sign matrix for classes lowDim to lowDim+1
	lowIndex := lowDim;
	highIndex := lowDim+1;
	
	orientationMatrix := RMatrixSpace(Rationals(), #representatives[highIndex], #representatives[lowIndex]) ! 0;	
	
	for i in [1..#representatives[highIndex]] do //Higher dim representatives
		if orientability[highIndex][i] then
			//Find orientable facets of the cell
			facets := facetsAsForms(representatives[highIndex][i]);
			
			cellBasis := formBasis(representatives[highIndex][i]); //Fixes an orientation of the vector space spanned by vertices of the cell
			
			for facet in facets do
				lowCellIndex := Position(cells[lowIndex], facet);
				facetIndex := equivalenceIndices[lowIndex][lowCellIndex];
				
				if orientability[lowIndex][facetIndex] then
					//Calculate relative orientation
					facetBasis := formBasis(facet); //Fixes an orientation of the vector space spanned by vertices of the facet
					
					for form in representatives[highIndex][i] do
						if Position(facet, form) eq 0 then
							newForm := form;
							break;
						end if;
					end for;
					
					M := MatrixRing(Rationals(), #cellBasis) ! linearCombinations(cellBasis, Append(facetBasis,newForm)); //Coercion necessary otherwise it segfaults
					relativeOrientation := Sign(Determinant(M));
					
					//Calculate orientation compatibility with representative
					repBasis := formBasis(representatives[lowIndex][facetIndex]);
					
					M := MatrixRing(Rationals(), #repBasis) ! linearCombinations(repBasis, [Dagger(equivalenceWitnesses[lowIndex][lowCellIndex]) * form * equivalenceWitnesses[lowIndex][lowCellIndex] : form in facetBasis]);
					orientationCompatibility := Sign(Determinant(M));
					
					//Contribution to cochain map
					orientationMatrix[i][facetIndex] +:= relativeOrientation * orientationCompatibility;
					//print relativeOrientation, orientationCompatibility, "from class", facetIndex, "and facet", lowCellIndex;
				end if;
			end for;
		end if;
	end for;
	
	return orientationMatrix;
end function;
