import "init.m" : matricesB, Dagger, n, overO, dieuDet, centraliserBasis;
import "symmetricSpace.m" : perpendicularForm, minimalVectors, positiveDefinite, getVectorsSizeRange, evaluateHermitian, perfectionRank, toHermitian, toHermitians, createGram;
import "polytope.m" : facetsAsForms;

function initialPerfectForm()
	//Start with the identity
	form := matricesB ! 1;
	
	min, minVecs := minimalVectors(form);
	rank := perfectionRank(form);
	
	//Now travel in orthgonal direction to add minimal vectors
	orthoForm := perpendicularForm(minVecs);
	
	while orthoForm ne 0 do //Still need to add another minimal vector
		searchSize := min; //We start by looking for vectors to add between the minimum and minimum+1
		
		found := false;
		while not found do
			addVecs := getVectorsSizeRange(form, searchSize, searchSize+1);
			
			for i in [1..#addVecs] do
				val := evaluateHermitian(orthoForm, addVecs[i]);
				
				if val ne 0 then
					newForm := form + (min - evaluateHermitian(form, addVecs[i]))/val * orthoForm; //newForm evaluates to min on minVecs and on addVecs[i]
					
					if positiveDefinite(newForm) then
						if  perfectionRank(newForm) gt rank then
							form := newForm;
							min, minVecs := minimalVectors(form);
							rank := perfectionRank(form);
							
							orthoForm := perpendicularForm(minVecs);
							found := true;
							break;
						end if;
					end if;
				end if;
			end for;
			
			searchSize +:= 1; //None of the above gave something positive definite? try again, except use vectors in a range 1 higher
		end while;
	end while;
	
	return form/min; //normalised to have minimum 1
end function;

function neighbour(perfectForm, facet) //finds the perfect form sharing facet with perfectForm
	orthoForm := perpendicularForm(facet);
	_, minVecs := minimalVectors(perfectForm);
	
	for i in [1..#minVecs] do
		if evaluateHermitian(orthoForm, minVecs[i]) lt 0 then //<orthoForm, minVecs[i]> >= 0 required
			orthoForm := -orthoForm;
		end if;
	end for;
	
	//First find an element l of S = {x in O^n : perfectForm + rho(x) * orthoForm positive definite; <orthoForm, x> < 0} (see Gunnells '99)
	searchSize := 1;
	found := false;
	
	while not found do
		addVecs := getVectorsSizeRange(perfectForm, searchSize, searchSize+1);
		
		for i in [1..#addVecs] do
			val := evaluateHermitian(orthoForm, addVecs[i]);
			if val lt 0 then
				rho := (1 - evaluateHermitian(perfectForm, addVecs[i]))/val;
				
				if positiveDefinite(perfectForm + rho * orthoForm) then
					l := perfectForm + rho * orthoForm;
					found := true;
					break;
				end if;
			end if;
		end for;
		
		searchSize +:= 1;
	end while;
	
	//Now find the finite set T_l := S intersect {<x, l> <= 1} - minimal vectors are those minimising rho on T_l
	T := getVectorsSizeRange(l, 0, 1); //Why are we not required to take negatives of vectors here?? something to run tests for correctness on
	
	min := Infinity();
	for i in [1..#T] do
		if evaluateHermitian(orthoForm, T[i]) ne 0 then
			rho := (1 - evaluateHermitian(perfectForm, T[i]))/evaluateHermitian(orthoForm, T[i]);
			
			if rho lt min then
				min := rho;
			end if;
		end if;
	end for;
	
	return perfectForm + min * orthoForm;
end function;

function clearDenoms(gram1, gram2)
	denoms := [];
	
	for i in [1..NumberOfRows(gram1)] do
		for j in [1..NumberOfColumns(gram1)] do
			Append(~denoms, Denominator(gram1[i][j]));
			Append(~denoms, Denominator(gram2[i][j]));
		end for;
	end for;
	
	denom := LCM(denoms);
	gram1 *:= denom;
	gram2 *:= denom;
	
	return gram1, gram2;
end function;

function equivalent(form1, form2)
	gram1 := createGram(form1);
	gram2 := createGram(form2);
	
	if Determinant(gram1) eq Determinant(gram2) then
		//Create twisted Gram matrices with integer coefficients
		gram1, gram2 := clearDenoms(gram1, gram2);
		
		twistedGrams1 := [gram1 * b : b in centraliserBasis];
		twistedGrams2 := [gram2 * b : b in centraliserBasis];
		
		for i in [1..#centraliserBasis] do
			twistedGrams1[i], twistedGrams2[i] := clearDenoms(twistedGrams1[i], twistedGrams2[i]);
		end for;
		
		twistedGrams1 := [MatrixRing(Integers(), 4*n) ! mat : mat in twistedGrams1];
		twistedGrams2 := [MatrixRing(Integers(), 4*n) ! mat : mat in twistedGrams2];
		
		return IsIsometric(twistedGrams1, twistedGrams2);
	else
		return false;
	end if;
end function;

function voronoiAlgorithm()
	perfectForms := [initialPerfectForm()];
	untested := perfectForms;
	
	print "Found initial perfect form";
	
	while #untested gt 0 do
		print #untested, "classes found but not yet tested;", #perfectForms, "total classes found";
		
		_, minVecs := minimalVectors(untested[1]);
		facets := facetsAsForms(toHermitians(minVecs));
		
		print #facets, "facets for current perfect form";
		for i in [1..#facets] do //Find all neighbouring forms
			print "Testing facet", i;
			form := neighbour(untested[1], facets[i]);
			
			//Check if neighbouring form is equivalent to an already known class
			new := true;
			for j in [1..#perfectForms] do
				if equivalent(form, perfectForms[j]) then
					new := false;
					break;
				end if;
			end for;
			
			if new then
				Append(~perfectForms, form);
				Append(~untested, form);
				print "New perfect form class found";
			end if;
		end for;
		
		Remove(~untested, 1);
	end while;
	
	return perfectForms;
end function;
