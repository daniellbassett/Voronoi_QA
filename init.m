//This files prepares all of the things about the quaternion algebra for use in other files
//To use: edit n for the dimension of matrices, a,b, for the definite rational quaternion algebra B, and O for the order of B

import "symmetricSpace.m" : innerProduct;

Q := Rationals(); //Base field

n := 2; //Number of variables
hermDim := n + 2*n*(n-1);

//Create quaternion algebra B/Q and order O
a := -2;
b := -5;
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
