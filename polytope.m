function toRationalSequence(A) //Converts a form over B into an element of V
	coords := [];
	
	for i in [1..NumberOfRows(A)] do
		for j in [1..NumberOfColumns(A)] do
			coords cat:= Coordinates(A[i][j]);
		end for;
	end for;
	
	return coords;
end function;

function toPolytope(S) //Takes in a sequence of forms (which will correspond to minimal vectors) and returns them as a polytope object
	vectors := [toRationalSequence(S[i]) : i in [1..#S]];
	
	return Polytope(vectors);
end function;

function facetsAsForms(S) //Takes a sequence of forms, returns a list of the facets of the polytope they form
	poly := toPolytope(S);
	facets := FacetIndices(poly);
	
	formList := [];
	for i in [1..#facets] do
		facet := IndexedSetToSequence(SetToIndexedSet(facets[i]));
		Append(~formList, [S[facet[j]] : j in [1..#facet]]);
	end for;
	
	return formList;
end function;
