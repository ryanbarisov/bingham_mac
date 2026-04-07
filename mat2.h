#ifndef _MAT2_H
#define _MAT2_H

#include <array>

// (0 2)
// (2 1)
//typedef double symmat2[3];//xx yy xy
typedef std::array<double,3> symmat2;
// (0 2)
// (1 3)
//typedef double mat2[4];//xx yx xy yy
typedef std::array<double,4> mat2;
// (0 3 6)
// (1 4 7)
// (2 5 8)
//typedef double mat3[9];//xx yx zx  xy yy zy xz yz zz
typedef std::array<double,9> mat3;

inline double det2(const symmat2& A) {return A[0]*A[1]-A[2]*A[2];}
inline double det3(const mat3& A) {return A[0]*(A[4]*A[8]-A[5]*A[7]) - A[3]*(A[1]*A[8]-A[2]*A[7]) + A[6]*(A[1]*A[5]-A[2]*A[4]);}
inline double tr(const symmat2& A) {return A[0]+A[1];}
// Frobenius norm of symmetric 2x2 matrix
inline double frobnorm(const symmat2& psi) {return sqrt(psi[0]*psi[0]+psi[1]*psi[1]+2*psi[2]*psi[2]);}

void eigs(const symmat2& A, double eigvals[2], mat2& eigvecs)
{
	static const bool scale = false;
	symmat2 A1 = A;
	double fnorm = frobnorm(A); if(scale && fnorm) { A1[0] /= fnorm; A1[1] /= fnorm; A1[2] /= fnorm; }
	double discr = tr(A1)*tr(A1) - 4*det2(A1), d1, d2;
	if(discr < -1.0e-15) throw "Negative discriminant, should not happen!";
	discr = std::max(0.0, discr);
	eigvals[0] = 0.5*(tr(A1) + sqrt(discr));
	eigvals[1] = 0.5*(tr(A1) - sqrt(discr));
	bool err = fabs(eigvals[0]-eigvals[1]) < 1.0e-12;
	if(!err)
	{
		eigvecs[0] = A1[2];
		eigvecs[1] = eigvals[0] - A1[0];
		eigvecs[2] = eigvals[1] - A1[1];
		eigvecs[3] = A1[2];
		d1 = 0;//sqrt(eigvecs[0]*eigvecs[0]+eigvecs[1]*eigvecs[1]);
		if(d1)
		{
			eigvecs[0] /= d1;
			eigvecs[1] /= d1;
		}
		else err = true;
		d2 = 0;//sqrt(eigvecs[2]*eigvecs[2]+eigvecs[3]*eigvecs[3]);
		if(d2)
		{
			eigvecs[2] /= d2;
			eigvecs[3] /= d2;
		}
		else err = true;
	}
	if(err)
	{
		eigvecs[0] = 1.0; eigvecs[1] = 0.0;
		eigvecs[2] = 0.0; eigvecs[3] = 1.0;
	}
	if(scale) {eigvals[0] *= fnorm; eigvals[1] *= fnorm;}
	if(fabs(eigvecs[0]*eigvecs[2] + eigvecs[1]*eigvecs[3]) > 1.0e-7)
		throw "non-orthogonal eigenvectors!";
}
inline void transpose(const mat2& A, mat2& At)
{
	At[0] = A[0];
	At[1] = A[2];
	At[2] = A[1];
	At[3] = A[3];
}
inline void matmul2ss(const symmat2& A, const symmat2& B, symmat2& AB)
{
	AB[0] = A[0]*B[0]+A[2]*B[2]; // xx
	AB[1] = A[2]*B[2]+A[1]*B[1]; // yy
	AB[2] = A[0]*B[2]+A[2]*B[1]; // xy
}
inline void matmul2ns(const mat2& A, const symmat2& B, mat2& AB)
{
	// ( A[0] A[2] ) ( B[0] B[2] )
	// ( A[1] A[3] ) ( B[2] B[1] )
	AB[0] = A[0]*B[0]+A[2]*B[2]; // xx
	AB[1] = A[1]*B[0]+A[3]*B[2]; // yx
	AB[2] = A[0]*B[2]+A[2]*B[1]; // xy
	AB[3] = A[1]*B[2]+A[3]*B[1]; // yy
}
inline void matmul2sn(const symmat2& A, const mat2& B, mat2& AB)
{
	// ( A[0] A[2] ) ( B[0] B[2] )
	// ( A[2] A[1] ) ( B[1] B[3] )
	AB[0] = A[0]*B[0]+A[2]*B[1]; // xx
	AB[1] = A[2]*B[0]+A[1]*B[1]; // yx
	AB[2] = A[0]*B[2]+A[2]*B[3]; // xy
	AB[3] = A[2]*B[2]+A[1]*B[3]; // yy
}
// B = R^T A R, B = B^T when A = A^T
inline void matmul3ns(const mat2& R, const symmat2& A, symmat2& RtAR)
{
	// ( R[0] R[2] ) -> ( r11 r12 ), ( A[0] A[2] )
	// ( R[1] R[3] ) -> ( r21 r22 ), ( A[2] A[1] )
	double r11 = R[0], r21 = R[1], r12 = R[2], r22 = R[3];
	double a11 = A[0], a21 = A[2], a12 = A[2], a22 = A[1];
	RtAR[0] = r11*r11*a11 + r11*r21*(a12+a21) + r21*r21*a22; // xx
	RtAR[1] = r12*r12*a11 + r12*r22*(a12+a21) + r22*r22*a22; // yy
	RtAR[2] = r11*r12*a11 + r11*r22*a12 + r12*r21*a21 + r21*r22*a22; // xy
}
// B = R^T A R
inline void matmul3nn(const mat2& R, const mat2& A, mat2& RtAR)
{
	// ( A[0] A[2] ) -> ( a11 a12 )
	// ( A[1] A[3] ) -> ( a21 a22 )
	double r11 = R[0], r21 = R[1], r12 = R[2], r22 = R[3];
	double a11 = A[0], a21 = A[1], a12 = A[2], a22 = A[3];
	RtAR[0] = r11*r11*a11 + r11*r21*(a12+a21) + r21*r21*a22; // xx
	RtAR[1] = r11*r12*a11 + r12*r21*a12 + r11*r22*a21 + r21*r22*a22; // yx
	RtAR[2] = r11*r12*a11 + r11*r22*a12 + r12*r21*a21 + r21*r22*a22; // xy
	RtAR[3] = r12*r12*a11 + r12*r22*(a12+a21) + r22*r22*a22; // yy
}
inline void solve3(mat3& A, const double b[3], double x[3])
{
	// Cramer method
	double v[3], d = det3(A);
	if(fabs(d) < 1.0e-12) throw "zero determinant";
	for(int i = 0; i < 3; ++i)
	{
		// overwrite column i of A with RHS vector
		for(int q = 0; q < 3; ++q)
		{
			v[q] = A[3*i+q];
			A[3*i+q] = b[q];
		}
		x[i] = det3(A)/d;
		// restore column i
		for(int q = 0; q < 3; ++q)
			A[3*i+q] = v[q];
	}
}


int test_mat2()
{
	symmat2 A1;
	A1[0] = 1e-2; A1[1] = -1e-2; A1[2] = 1e-5;
	mat2 eigvecs;
	double eigvals[2];
	eigs(A1, eigvals, eigvecs);
	std::cout << "eigvals: b1 " << eigvals[0] << " b2 " << eigvals[1] << std::endl;
	std::cout << "eigvecs: v1 " << eigvecs[0] << " " << eigvecs[1] << std::endl;
	std::cout << "eigvecs: v2 " << eigvecs[2] << " " << eigvecs[3] << std::endl;
	mat2 AV;
	matmul2sn(A1, eigvecs, AV);
	std::cout << "matvec: A*v1 " << AV[0] << " " << AV[1]
		<< " A*v2 " << AV[2] << " " << AV[3] << std::endl;
	std::cout << "matvec: A*v1-b1*v1 " << (AV[0]-eigvals[0]*eigvecs[0]) << " " << (AV[1]-eigvals[0]*eigvecs[1]) << std::endl;
	std::cout << "matvec: A*v2-b2*v2 " << (AV[2]-eigvals[1]*eigvecs[2]) << " " << (AV[3]-eigvals[1]*eigvecs[3]) << std::endl;

	mat2 A, AB, C, G, AGAt, H, At;
	symmat2 B, D, E, DB, ABAt, F;
	A[0] = 1; A[1] = 3; A[2] = 2; A[3] = 4;
	B[0] = 5; B[1] = 8; B[2] = 6;
	G[0] = 5; G[1] = 6; G[2] = 7; G[3] = 8;
	C[0] = 17; C[1] = 39; C[2] = 22; C[3] = 50;
	F[0] = 61; F[1] = 317; F[2] = 139;
	D[0] = 1; D[1] = 1; D[2] = 2;
	E[0] = 17; E[1] = 20; E[2] = 22;
	H[0] = 63; H[1] = 145; H[2] = 143; H[3] = 329;
	matmul2ns(A,B,AB);
	std::cout << "A*B - AB = " << (AB[0]-C[0]) << " " << (AB[1]-C[1]) << " " << (AB[2]-C[2]) << " " << (AB[3]-C[3]) << std::endl;
	matmul2ss(D,B,DB);
	std::cout << "D*B - DB = " << (DB[0]-E[0]) << " " << (DB[1]-E[1]) << " " << (DB[2]-E[2]) << std::endl;

	transpose(A,At);
	matmul3nn(At,G,AGAt);
	std::cout << "A*G*At - AGAt = " << (AGAt[0]-H[0]) << " " << (AGAt[1]-H[1]) << " " << (AGAt[2]-H[2]) << " " << (AGAt[3]-H[3]) << std::endl;

	matmul3ns(At,B,ABAt);
	std::cout << "A*B*At - ABAt = " << (ABAt[0]-F[0]) << " " << (ABAt[1]-F[1]) << " " << (ABAt[2]-F[2]) << std::endl;

	std::cout << "det2 is ok? " << (fabs(det2(B)-4) < 1.0e-12) << std::endl;

	mat3 A3 = {1,0,0,2,5,6,3,6,8};
	std::cout << "det3 is ok? " << (fabs(det3(A3)-4) < 1.0e-12) << std::endl;

	double rhs[3] = {14,28,36}, x[3] = {1,2,3}, sol[3];
	solve3(A3, rhs, sol);
	std::cout << "solve3: sol " << sol[0] << " " << sol[1] << " " << sol[2] << "; expected " << x[0] << " " << x[1] << " " << x[2] << std::endl;

	return 0;
}

#endif //_MAT2_H
