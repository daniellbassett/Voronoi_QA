//This files prepares all of the things about the quaternion algebra for use in other files
//To use: edit n for the dimension of matrices, a,b, for the definite rational quaternion algebra B, and O for the order of B

import "symmetricSpace.m" : innerProduct;

Q := Rationals(); //Base field

n := 2; //Number of variables
hermDim := n + 2*n*(n-1);

//Create quaternion algebra B/Q and order O
a := -1;
b := -1;
B<i,j,k> := QuaternionAlgebra<Q|a,b>;
O := MaximalOrder(B); 

//Basis for order
orderBasis := Basis(O);
matricesB := MatrixRing(B, n);

//Create basis for O^n from a basis for O, as a list of vectors
latticeBasis := [RMatrixSpace(B, n, 1) ! 0 : i in [1..4*n]];
for i in [1..n] do
	for j in [1..4] do
		latticeBasis[j+4*(i-1)][i][1] := orderBasis[j];
	end for;
end for;

//Involution on matrices over B
function Dagger(A)
	ADagger := RMatrixSpace(B, NumberOfColumns(A), NumberOfRows(A)) ! 0;
	
	for i in [1..NumberOfRows(A)] do
		for j in [1..NumberOfColumns(A)] do
			ADagger[j][i] := Conjugate(A[i][j]);
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
		for k in [1..4] do //orderBasis has length 4 since quaternion algebra
			hermBasis[count][i][j] := orderBasis[k];
			hermBasis[count][j][i] := Conjugate(orderBasis[k]);
			
			count +:= 1;
		end for;
	end for;
end for;

//Centraliser of the embedding as rational matrices, for use in equivalence testing and automorphism group calculations (Coulangeon et. al. Lemma 7.2)
orderBasisCoordinates := Transpose(MatrixRing(Q,4) ! [Coordinates(orderBasis[i]) : i in [1..4]]);

function leftRegularRep(x) //x acting on the right on B as a Q-vector space
	coords := [];
	
	for w in orderBasis do
		Append(~coords, RMatrixSpace(Q, 1, 4) ! Coordinates(x*w));
	end for;
	
	mat := Transpose(MatrixRing(Rationals(), 4) ! coords);
	return orderBasisCoordinates^-1 * mat;
end function;

function rationalMatrixEmbedding(A)
	mat := MatrixRing(Rationals(), 4*n) ! 0;
	for i in [1..n] do
		for j in [1..n] do
			elementEmbedding := leftRegularRep(A[i][j]);
			
			for k in [1..4] do
				for l in [1..4] do
					mat[4*(i-1)+k][4*(j-1)+l] := elementEmbedding[k][l];
				end for;
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

embeddedAlgebra := sub<MatrixAlgebra(Rationals(), 4*n) | embeddedMatrices>;
embeddedCentraliser := Centraliser(MatrixAlgebra(Rationals(), 4*n), embeddedAlgebra);
centraliserBasis := Basis(embeddedCentraliser);

function matrixRationalToQuaternion(A)
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

//Dieudonne determinant of 2x2 matrices
function dieuDet(A)
	return Norm(A[1][1] * A[2][2]) + Norm(A[1][2]*A[2][1]) - Trace(A[1][1] * Conjugate(A[2][1]) * A[2][2] * Conjugate(A[1][2]));
end function;
