//The symmetric space of Hermitian nxn matrices over a quaternion algebra
import "init.m" : B, n, latticeBasis, hermBasis, Dagger;


function innerProduct(A,B)
	return Trace(Trace(A * B)); //Inner trace is matrix trace, outer trace is quaternion algebra reduced trace; this composition is the reduced trace on Mat_n(H)/R
end function;

//Embeds O^n vectors in Herm_n
function toHermitians(S)
	forms := [];
	for i in [1..#S] do
		form := S[i] * Dagger(S[i]);
		new := true;
		
		for j in [1..#forms] do
			if form eq forms[j] then
				new := false;
				break;
			end if;
		end for;
		
		if new then
			Append(~forms, form); //i think ~ is a pointer thing
		end if;
	end for;
	
	return forms;
end function;


//-----Finding minimal vectors and minima-----

//Evaluating forms
function evaluateHermitian(A, v) //Evaluates Hermitian form A on vector v
	return innerProduct(A, v*Dagger(v));
end function;

function evaluateBilinear(A, v, w) //Evaluates the corresponding real bilinear form on vectors v,w
	return Trace(Trace(A*(v*Dagger(w) + w*Dagger(v))))/2;
end function;

//Gram matrix
function createGram(A)
	Gram := MatrixRing(Rationals(), #latticeBasis) ! 0; //Rationals() instead of the number field version allows use of lattice reduction routines, since it is considered as a subfield of the reals
	
	for i in [1..#latticeBasis] do
		Gram[i][i] := evaluateHermitian(A, latticeBasis[i]);
		
		for j in [i+1..#latticeBasis] do
			Gram[i][j] := evaluateBilinear(A, latticeBasis[i], latticeBasis[j]);
			Gram[j][i] := Gram[i][j];
		end for;
	end for;
	
	return Gram;
end function;

function positiveDefinite(A)
	Gram := createGram(A);
	return IsPositiveDefinite(Gram);
end function;


//Minima and minimal vectors
function minimalVectors(A)
	Gram := createGram(A);
	L := LatticeWithGram(Gram);
	
	coefficients := ShortestVectors(L); //gets coefficients of the lattice for the minimal vectors; also assigns L`minimum
	
	minVecs := [RMatrixSpace(B, n, 1) ! 0 : i in [1..2*#coefficients]]; //ShortestVectors only returns minimal vectors up to sign; we add negatives back in
	for i in [1..#coefficients] do
		for j in [1..#latticeBasis] do
			minVecs[i] +:=  coefficients[i][j] * latticeBasis[j];
			minVecs[i+#coefficients] +:= -coefficients[i][j] * latticeBasis[j]; 
		end for;
	end for;
	
	return L`Minimum, minVecs;
end function;

function getVectorsSizeRange(A, lower, upper)
	Gram := createGram(A);
	L := LatticeWithGram(Gram);
	
	if lower eq 0 then
		coefficients := ShortVectors(L, upper);
	else
		coefficients := ShortVectors(L, lower, upper);
	end if;
	
	shortVecs := [RMatrixSpace(B, n, 1) ! 0 : i in [1..#coefficients]];
	for i in [1..#coefficients] do
		for j in [1..#latticeBasis] do
			shortVecs[i] +:= coefficients[i][1][j] * latticeBasis[j];
		end for;
	end for;
	
	return shortVecs;
end function;

//Perpendicular forms
function perpendicularForm(S) //A set S of vectors in O^n
	innerProductMatrix := RMatrixSpace(B, #hermBasis, #S) ! 0;

	for i in [1..#hermBasis] do
		for j in [1..#S] do
			innerProductMatrix[i][j] := innerProduct(hermBasis[i], S[j] * Dagger(S[j])); //The inner product with the corresponding forms to the vectors
		end for;
	end for;
	
	kernelCoordinates := KernelMatrix(innerProductMatrix);
	
	orthogonalForm := MatrixRing(B, n) ! 0;
	
	if NumberOfRows(kernelCoordinates) gt 0 then //At least one orthogonal form found
		for i in [1..#hermBasis] do
			orthogonalForm +:= kernelCoordinates[1][i] * hermBasis[i];
		end for;
	end if;
	
	return orthogonalForm;
end function;

function perfectionRank(form)
	_, minVecs := minimalVectors(form);
	S := toHermitians(minVecs);
	innerProductMatrix := RMatrixSpace(B, #hermBasis, #S) ! 0;

	for i in [1..#hermBasis] do
		for j in [1..#S] do
			innerProductMatrix[i][j] := innerProduct(hermBasis[i], S[j] * Dagger(S[j])); //The inner product with the corresponding forms to the vectors
		end for;
	end for;
	
	kernelCoordinates := KernelMatrix(innerProductMatrix);
	return #hermBasis - NumberOfRows(kernelCoordinates);
end function;
