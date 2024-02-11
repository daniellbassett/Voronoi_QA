import "init.m" : matricesB, elementToRational;

//Creating polytopes from forms
function toRationalSequence(A) //Converts a form over B into an element of V
	coords := [];
	
	for i in [1..NumberOfRows(A)] do
		for j in [1..NumberOfColumns(A)] do
			coords cat:= elementToRational(A[i][j]);
		end for;
	end for;
	
	return coords;
end function;

function toPolytope(S) //Takes in a sequence of forms (which will correspond to minimal vectors) and returns them as a polytope object
	vectors := [toRationalSequence(S[i]) : i in [1..#S]];
	
	return Polytope(vectors);
end function;

//Finding facets
function facetsAsForms(S) //Takes a sequence of forms, returns a list of the facets of the polytope they form
	poly := toPolytope(S);
	
	redundant := true;
	while redundant do
		redundant := false;
		for i in [1..#S] do
			otherForms := Remove(S, i);
			poly2 := toPolytope(otherForms);
			
			if poly eq poly2 then
				redundant := true;
				Remove(~S, i);
				break;
			end if;
		end for;
	end while;
	
	poly := toPolytope(S);
	
	facets := FacetIndices(poly);
	
	formList := [];
	for i in [1..#facets] do
		facet := IndexedSetToSequence(SetToIndexedSet(facets[i]));
		Append(~formList, [S[facet[j]] : j in [1..#facet]]);
	end for;
	
	return formList;
end function;

//Barycentre of a facet defined by forms
function barycentre(S)
	sum := matricesB ! 0;
	for form in S do
		sum +:= form;
	end for;
	
	return sum/#S;
end function;
