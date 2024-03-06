d := 427;

R<t> := PolynomialRing(Integers());
K<theta> := NumberField(t^2+d); //Base field
sigma := Automorphisms(K)[2];

n := 2; //Number of variables
hermDim := n + n*(n-1); //Each diagonal entry adds 1, each non-diagonal adds 2

//Order and basis
O := MaximalOrder(K); 
orderBasis := Basis(O);

matricesB := MatrixRing(K, n);

//Create basis for O^n from a basis for O, as a list of vectors
latticeBasis := [RMatrixSpace(K, n, 1) ! 0 : i in [1..#orderBasis*n]];
for i in [1..n] do
	for j in [1..#orderBasis] do
		latticeBasis[j+#orderBasis*(i-1)][i][1] := orderBasis[j];
	end for;
end for;

//Involution on matrices over B
function Dagger(A)
	ADagger := RMatrixSpace(K, NumberOfColumns(A), NumberOfRows(A)) ! 0;
	
	for i in [1..NumberOfRows(A)] do
		for j in [1..NumberOfColumns(A)] do
			ADagger[j][i] := sigma(A[i][j]);
		end for;
	end for;
	
	return ADagger;
end function;

//Create basis for Herm(B^n) from a basis for O, as a list of matrices
hermBasis := [matricesB ! 0 : i in [1..hermDim]];
count := 1;
for i in [1..n] do
	hermBasis[count][i][i] := 1;
	count +:= 1;
	
	for j in [i+1..n] do
		for k in [1..#orderBasis] do
			hermBasis[count][i][j] := orderBasis[k];
			hermBasis[count][j][i] := sigma(orderBasis[k]);
			
			count +:= 1;
		end for;
	end for;
end for;

//Centraliser of the embedding as rational matrices, for use in equivalence testing and automorphism group calculations (Coulangeon et. al. Lemma 7.2)
/*function elementToRational(x)
	coords := Coordinates(x);
	coordsRat := [];
	for y in coords do
		coordsRat cat:= Eltseq(y);
	end for;
	
	return coordsRat;
end function;*/

orderBasisCoordinates := Transpose(MatrixRing(Rationals(),#orderBasis) ! [Eltseq(K ! orderBasis[i]) : i in [1..#orderBasis]]);

function leftRegularRep(x) //x acting on the right on B as a Q-vector space
	coords := [];
	
	for w in orderBasis do
		Append(~coords, RMatrixSpace(Rationals(), 1, #orderBasis) ! Eltseq(x*w));
	end for;
	
	mat := Transpose(MatrixRing(Rationals(), #orderBasis) ! coords);
	return orderBasisCoordinates^-1 * mat;
end function;

function rationalMatrixEmbedding(A)
	mat := MatrixRing(Rationals(), #orderBasis*n) ! 0;
	for i in [1..n] do
		for j in [1..n] do
			elementEmbedding := leftRegularRep(A[i][j]);
			
			for k in [1..#orderBasis] do
				for l in [1..#orderBasis] do
					mat[#orderBasis*(i-1)+k][#orderBasis*(j-1)+l] := elementEmbedding[k][l];
				end for;
			end for;
		end for;
	end for;
	
	return mat;
end function;

function matrixRationalToQuaternion(A)
	mat := matricesB ! 0;
	
	for i in [1..n] do
		for j in [1..n] do
			for k in [1..#orderBasis] do
				mat[i][j] +:= orderBasis[k] * A[#orderBasis*(i-1)+k][#orderBasis*(j-1)+1]; //Regular representation was taken in basis orderBasis
			end for;
		end for;
	end for;
	
	return mat;
end function;

embeddedMatrices := [];
for i in [1..n] do
	for j in [1..n] do
		for w in orderBasis do
			mat := matricesB ! 0;
			mat[i][j] := w;
			
			Append(~embeddedMatrices, rationalMatrixEmbedding(mat));
		end for;
	end for;
end for;

embeddedAlgebra := sub<MatrixAlgebra(Rationals(), #orderBasis*n) | embeddedMatrices>;
embeddedCentraliser := Centraliser(MatrixAlgebra(Rationals(), #orderBasis*n), embeddedAlgebra);
centraliserBasis := Basis(embeddedCentraliser);

//Checks if matrix A is defined over the order O
function overO(A)
	for i in [1..NumberOfRows(A)] do
		for j in [1..NumberOfColumns(A)] do
			if not A[i][j] in O then
				return false;
			end if;
		end for;
	end for;
	
	return true;
end function;

function matrixRationalToImaginary(A)
	mat := matricesB ! 0;

	for i in [1..n] do
		for j in [1..n] do
			for k in [1..#orderBasis] do
				mat[i][j] +:= orderBasis[k] * A[#orderBasis*(i-1)+k][#orderBasis*(j-1)+1];
			end for;
		end for;
	end for;
	
	return mat;
end function;
/*
//Dieudonne determinant of 2x2 matrices
function dieuDet(A)
	return Norm(A[1][1] * A[2][2]) + Norm(A[1][2]*A[2][1]) - Trace(A[1][1] * Conjugate(A[2][1]) * A[2][2] * Conjugate(A[1][2]));
end function;*/
