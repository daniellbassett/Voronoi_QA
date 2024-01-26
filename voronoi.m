import "init.m" : matricesB, Dagger, n, O;
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
	T := getVectorsSizeRange(l, 0, 1);
	
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

function equivalent(form1, form2)
	//Currently hard-coded for n = 2.
	_, minVecs1 := minimalVectors(form1);
	_, minVecs2 := minimalVectors(form2);
	
	if #minVecs1 eq #minVecs2 then
		if Determinant(createGram(form1)) eq Determinant(createGram(form2)) then
			//First find linearly independent minimal vectors of form1
			found1 := false;
			for i in [1..#minVecs1] do
				for j in [i+1..#minVecs1] do
					if toHermitian(minVecs1[i]) ne toHermitian(minVecs1[j]) then
						found1 := true;
						
						i1 := i;
						j1 := j; 
						break;
					end if;
				end for;
				
				if found1 then
					break;
				end if;
			end for;
			
			mat1 := Transpose(matricesB ! [Transpose(minVecs1[i1]), Transpose(minVecs1[j1])]);
			mat1Inv := mat1^-1;
			
			//Now find linear map taking them to each pair of minimal vectors of form2; if form1 and form2 are equivalent, it must be under a map of this form
			for i in [1..#minVecs2] do
				for j in [1..#minVecs2] do //Need to go from 1 instead of i+1 since the linear map may not preserve the ordering of the vertices
					//Check the minimal vectors of form 2 are also independent
					if toHermitian(minVecs2[i]) ne toHermitian(minVecs2[j]) then
						//Create linear map
						mat2 := Transpose(matricesB ! [Transpose(minVecs2[i]), Transpose(minVecs2[j])]);
						mat := mat2 * mat1Inv;
						
						det := Determinant(mat);
						if det ne 0 then
							//Is it defined over O
							overO := true;
							
							for k in [1..n] do
								for l in [1..n] do
									if not mat[k][l] in O then
										overO := false;
										break;
									end if;
								end for;
								
								if not overO then
									break;
								end if;
							end for;
							
							if overO then
								//Is the inverse defined over O
								matInv := mat^-1;
								
								for k in [1..n] do
									for l in [1..n] do
										if not matInv[k][l] in O then
											overO := false;
											break;
										end if;
									end for;
									
									if not overO then
										break;
									end if;
								end for;
								
								if overO then
									//Does it take form2 to form1 as desired
									if form1 eq Dagger(mat)*form2*mat then
										return true;
									end if;
								end if;
							end if;
						end if;
					end if;
				end for;
			end for;
		end if;
	end if;
	
	return false;
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
