//[[Rcpp::depends(RcppEigen)]]

#include "tensorTn.h"

//***--------------------updateS in T4----------------------------------**
MatrixXd updateT4S(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D)
{
	int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, r4 = opts.r4;
	int d = r1*r2*r3*r4, k,k1,j,j1;
	MatrixXd ztilde = Z * kroneckerProduct(C,kroneckerProduct(B, A)), vectorS;
	VectorXd U;
	U.setZero(d);
	MatrixXd  V = MatrixXd::Constant(d, d, 0);

	for (k = 0; k < r1*r2*r3; k++) {
		for (j = 0; j < r4; j++) {
			U[k*r4 + j] = ztilde.col(k).transpose()*Y*D.col(j);
			for (k1 = 0; k1 < r1*r2*r3; k1++) {
				for (j1 = 0; j1 < r4; j1++) {
					V(k*r4 + j, k1*r4 + j1) = kroneckerProduct(
						ztilde.col(k1).array()*ztilde.col(k).array(),
						(D.col(j1).array()*D.col(j).array()).transpose()).sum();
				}
			}
		}
	}
	vectorS = V.colPivHouseholderQr().solve(U);
	vectorS.resize(r4, r1*r2*r3);
	return vectorS;
}
//----------------------------------------------------------------**
//***--------------------updateD in T4----------------------------------**
MatrixXd updateT4D(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D)
{
	int q = opts.q,j,kp;
	MatrixXd ztilde = Z * kroneckerProduct(C,kroneckerProduct(B, A));
	MatrixXd StZ = ztilde * S.transpose();
	MatrixXd Cnew = MatrixXd::Constant(q, opts.r4, 0);
	HouseholderQR<MatrixXd> qr;
	qr.compute(StZ);
	MatrixXd R = qr.matrixQR().triangularView<Upper>();
	MatrixXd Q = qr.householderQ();

	
	kp = StZ.cols();
	MatrixXd temp, IDEN = MatrixXd::Identity(kp, kp);
	if (pow(condition_numberQRSym(R),2) > 1e10){
	  temp = tRbyR(R) + (IDEN.array()*1e-4).matrix();
	  for (j = 0; j < q; j++) Cnew.row(j) = (temp.colPivHouseholderQr().solve(StZ.transpose()*Y.col(j))).transpose();
	}
	else
	  for (j = 0; j < q; j++) Cnew.row(j) = (QbyR(Q.transpose(), UpTriangularInv(R), 0)*Y.col(j)).transpose();
	
	return Cnew;
}
//----------------------------------------------------------------**
//***--------------------updateA in T4----------------------------------**
MatrixXd updateT4A(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D)
{
	int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, K = opts.K, G = opts.G, p = opts.p;
	int d = r1 * p, t1,t2,t3,t4,j;
	MatrixXd W = D * S, Wt1, Zt2, Wt3, Zt4, zbw1, zbw2, vectorA, Gamma = kroneckerProduct(C,B);
	VectorXd tU;
	tU.setZero(d);
	MatrixXd tV = MatrixXd::Constant(d, d, 0);

	for (t2 = 0; t2<p; t2++) {
		Zt2 = Z.col(t2);
		for (j = 1; j<K*G; j++)	Zt2 = cbind_rcpp(Zt2, Z.col(j*p + t2));
		for (t1 = 0; t1<r1; t1++) {
			Wt1 = W.col(t1);
			for (j = 1; j<r2*r3; j++) Wt1 = cbind_rcpp(Wt1, W.col(j*r1 + t1));
			zbw1 = Zt2 * Gamma *(Wt1.transpose());
			tU[t2*r1 + t1] = (Y.array()*zbw1.array()).sum();
			for (t4 = 0; t4<p; t4++) {
				Zt4 = Z.col(t4);
				for (j = 1; j<K*G; j++) Zt4 = cbind_rcpp(Zt4, Z.col(j*p + t4));
				for (t3 = 0; t3<r1; t3++) {
					Wt3 = W.col(t3);
					for (j = 1; j<r2*r3; j++)	Wt3 = cbind_rcpp(Wt3, W.col(j*r1 + t3));
					zbw2 = Zt4 * Gamma *(Wt3.transpose());
					tV(t2*r1 + t1, t4*r1 + t3) = (zbw1.array()*zbw2.array()).sum();
				}
			}
		}
	}
	vectorA = tV.colPivHouseholderQr().solve(tU);
	vectorA.resize(r1, p);
	return vectorA.transpose();
}
//----------------------------------------------------------------**
//***--------------------updateB in T4----------------------------------**
MatrixXd updateT4B(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D)
{
    int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, K = opts.K, G = opts.G, q = opts.q, n = opts.n, p = opts.p; // p = A.rows();
	int d = r2 * K, t1, t2, t3, t4, j;
	MatrixXd W = D * S, Wt1, Zt2, Wt3, Zt4, zaw1, zaw2, vectorB, Gamma = kroneckerProduct(C,A);
	VectorXd tU;
	tU.setZero(d);
    MatrixXd tV = MatrixXd::Constant(d, d, 0);
    for (t2 = 0; t2<K; t2++) {
        Zt2 = Z.block(0, t2*p, n, p);
        for (j = 1; j<G; j++)    Zt2 = cbind_rcpp(Zt2, Z.block(0, j*p*K + t2*p, n, p));
        for (t1 = 0; t1<r2; t1++) {
            Wt1 = W.block(0, t1*r1, q, r1);
            for (j = 1; j<r3; j++) Wt1 = cbind_rcpp(Wt1, W.block(0, j*r1*r2 + t1*r1, q, r1));
            zaw1 = Zt2 * Gamma * (Wt1.transpose());
            tU[t2*r2 + t1] = (Y.array()*zaw1.array()).sum();
            for (t4 = 0; t4<K; t4++) {
                Zt4 = Z.block(0, t4*p, n, p);
                for (j = 1; j<G; j++)    Zt4 = cbind_rcpp(Zt4, Z.block(0, j*p*K + t4*p, n, p));
                for (t3 = 0; t3<r2; t3++) {
                    Wt3 = W.block(0, t3*r1, q, r1);
                    for (j = 1; j<r3; j++) Wt3 = cbind_rcpp(Wt3, W.block(0, j*r1*r2 + t3*r1, q, r1));
                    zaw2 = Zt4 * Gamma * (Wt3.transpose());
                    tV(t2*r2 + t1, t4*r2 + t3) = (zaw1.array()*zaw2.array()).sum();
                }
            }
        }
    }
	vectorB = solveEquationQR(tV, tU);
	vectorB.resize(r2, K);
	return vectorB.transpose();
}
//----------------------------------------------------------------**
//***--------------------updateC in T4----------------------------------**
MatrixXd updateT4C(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D)
{
    int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, K = opts.K, G = opts.G, q = opts.q, n = opts.n, p = opts.p;
	int d = r3 * G, t1, t2, t3, t4, j;
	MatrixXd W = D * S, Wt1, Zt2, Wt3, Zt4, zaw1, zaw2, vectorB, Gamma = kroneckerProduct(B,A);
	VectorXd tU;
	tU.setZero(d);
	MatrixXd tV = MatrixXd::Constant(d, d, 0);
    
    for (t2 = 0; t2<G; t2++) {
        Zt2 = Z.block(0, t2*p*K, n, p*K);
        for (t1 = 0; t1<r3; t1++) {
            Wt1 = W.block(0, t1*r1*r2, q, r1*r2);
            zaw1 = Zt2 * Gamma *(Wt1.transpose());
            tU[t2*r3 + t1] = (Y.array()*zaw1.array()).sum();
            for (t4 = 0; t4<G; t4++) {
                Zt4 = Z.block(0, t4*p*K, n, p*K);
                for (t3 = 0; t3<r3; t3++) {
                    Wt3 = W.block(0, t3*r1*r2, q, r1*r2);
                    zaw2 = Zt4 * Gamma *(Wt3.transpose());
                    tV(t2*r3 + t1, t4*r3 + t3) = (zaw1.array()*zaw2.array()).sum();
                }
            }
        }
    }
//    vectorB = solveEquationQR(tV, tU);
    vectorB = tV.colPivHouseholderQr().solve(tU);
	vectorB.resize(r3, G);
	return vectorB.transpose();
}
//----------------------------------------------------------------**
//***--------------------T4 Estimation without penalty---------------**
// [[Rcpp::export]]
List EstimationT4(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, List optsList)
{
	/*
	Input:
	Y is n*q matrix
	Z is n*(K*p*G) matrix
	A is p*r1 matrix
	B is K*r2 matrix
	C is G*r3 matrix
	D is q*r4 matrix
	S is r4*(r1*r2*r3) matrix
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	gamma is another tuning parameter for MCP and SCAD
	penalty= 1(lasso),2(mcp), and 3(scad)
	dfmax is the preset maximum digrees of freedom
	threshold is the error to control convergence for outer iteration
	eps is the error to control convergence for inner iteration
	max_step is the max step to control convergence for outer iteration
	max_iter is the max step to control convergence for inner iteration
	is_setlam is logical, 1 for set lambda by data; 0 for given lambda
	setlam is a vector saving preseting lam_max,lam_min,alpha,nlam
	nlam is the number of tuning parameters

	Output:
	Dnew is a estimator of D(4) = D%*%S%*%kronecker(t(C),kronecker(t(B),t(A)))
	*/

	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);
	opts.max_step1 = as<int>(optsList["max_step1"]);
	opts.r1 = as<int>(optsList["r1"]);
	opts.r2 = as<int>(optsList["r2"]);
	opts.r3 = as<int>(optsList["r3"]);
	opts.r4 = as<int>(optsList["r4"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.G = as<int>(optsList["G"]);
	opts.degr = as<int>(optsList["degr"]);
	opts.K = as<int>(optsList["K"]);
//    opts.n = Y.cols();
    opts.n = Y.rows();
	
	int j, max_step=opts.max_step;
	double  likhd0 = pow(10, 6), likhd1, eps = opts.eps;
	MatrixXd Snew, Anew, Bnew, Cnew, Dnew, Dn;
	VectorXi convergence1 = VectorXi::Constant(5, 1);
//    VectorXi activeA = VectorXi::Constant(p, 1);
//    MatrixXi betapath = MatrixXi::Constant(p, nlam, 0);
	
	int step = 0;
	while(step < max_step){		
		step = step + 1;
		Snew = updateT4S(Y, Z, A, B, C, D);
		Dn = D * Snew*kroneckerProduct(C,kroneckerProduct(B, A)).transpose();
		likhd1 = (Y - Z * Dn.transpose()).squaredNorm();
		Rcout << "step= " << step << "        S likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
		if (likhd1<likhd0) {
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[0]=0;	
		Dnew = updateT4D(Y, Z, S, A, B, C, D);
		Dn = Dnew * S*kroneckerProduct(C,kroneckerProduct(B, A)).transpose();
		likhd1 = (Y - Z * Dn.transpose()).squaredNorm();	
		Rcout << "step= " << step << "        D likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
		if (likhd1<likhd0) {
			D = Dnew;
			likhd0 = likhd1;
		}
		else convergence1[1]=0;
		Cnew = updateT4C(Y, Z, S, A, B, C, D);
		Dn = D * S* kroneckerProduct(Cnew,kroneckerProduct(B, A)).transpose();
		likhd1 = (Y - Z * Dn.transpose()).squaredNorm();	
		Rcout << "step= " << step << "        C likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
		if (likhd1<likhd0) {
			C = Cnew;
			likhd0 = likhd1;
		}
		else convergence1[2]=0;		
		Bnew = updateT4B(Y, Z, S, A, B, C, D);
		Dn = D * S* kroneckerProduct(C,kroneckerProduct(Bnew, A)).transpose();
		likhd1 = (Y - Z * Dn.transpose()).squaredNorm();
        Rcout << "step= " << step << "        B likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;		
		if (likhd1<likhd0) {
			B = Bnew;
			likhd0 = likhd1;		
		}
		else convergence1[3]=0;	
		Anew = updateT4A(Y, Z, S, A, B, C, D);
		Dn = D * S* kroneckerProduct(C,kroneckerProduct(B, Anew)).transpose();
		likhd1 = (Y - Z * Dn.transpose()).squaredNorm();	
        Rcout << "step= " << step << "        A likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;		
        //Rcout << "esp= " << opts.eps << "        fabs(likhd0-likhd1)/likhd0 = "<< fabs(likhd0-likhd1)/likhd0 << std::endl;
		if (likhd1<likhd0) {
			A = Anew;
			if(fabs(likhd0-likhd1)/fabs(likhd0+1) < eps) break;
			likhd0 = likhd1;			
		}
		else convergence1[4]=0;	
		//Rcout << "opts.max_step = " << opts.max_step << "        convergence1.sum()= " << convergence1.sum()<< std::endl;
		if(convergence1.sum()==0) break;		
		else for(j=0;j<5;j++) convergence1[j] = 1;
		//Rcout << "step= " << step << "        likhd1= " << likhd1<< std::endl;
//        if (step==5) stop("stop for test!");
	}
	return List::create(Named("likhd") = likhd1, Named("Snew") = S, Named("Anew") = A, 
	                    Named("Bnew") = B, Named("Cnew") = C, Named("Dnew") = D);
}
//----------------------------------------------------------------**
//***--------------------setup tuning parameters------------------**
// [[Rcpp::export]]
VectorXd setuplambdaT4(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, int nlam, VectorXd setlam)
{
	int n = Y.rows(), q = Y.cols(), G = C.rows(), p = A.rows()/G, K = B.rows(), j, jj, KG=K*G;
	int r1 = A.cols(), r2 = B.cols(), r3 = C.cols(), r4 = D.cols();
	double lam_max=setlam[0], lam_min=setlam[1], alpha=setlam[2];
	VectorXd lambda, lambda1, tmp, tmp1;
	MatrixXd S1, dcbs, V, V_1, Y1 = Y, Gammaj, Gamma_sqrt, svdu, svdd, U, L;
	Y1.resize(n*q, 1);
    VectorXi dims = VectorXi::Constant(4, 0);
    dims[0] = r1; dims[1] = r2; dims[2] = r3; dims[3] = r4;
	S1 = TransferModalUnfoldingsT(S, 4, 1, dims);
	dcbs = kroneckerProduct(D,kroneckerProduct(C, B))*(S1.transpose());
	VectorXi id = SEQ(1, p*KG, p);
	tmp = VectorXd::Constant(p, 0);
	for (j = 0; j < p; j++)
	{
		V = submatrix_col(Z, id.array() + j)*(dcbs.block(0, 0, KG, r1));
		for (jj = 1; jj < q; jj++) {
			V_1 = submatrix_col(Z, id.array() + j)*(dcbs.block(jj*KG, 0, KG, r1));
			V = rbind_rcpp(V, V_1);
		}
		Gammaj = ((V.transpose()*V).array() / n).matrix();		
		
		JacobiSVD<MatrixXd> svd(Gammaj, ComputeThinU | ComputeThinV);
		svdu = svd.matrixU();
		svdd = (svd.singularValues()).asDiagonal();
		Gamma_sqrt = svdu * ((1 / (svdd.diagonal().array().sqrt())).matrix().asDiagonal())*(svdu.transpose());
		tmp1 = Y1.transpose()* V * Gamma_sqrt;
		tmp[j]=tmp1.array().abs().sum();
	}

	double max_tmp;
	//max_tmp = ((tmp.array()).maxCoeff())*log(p)/sqrt(n*q);
	max_tmp = (tmp.array()).maxCoeff()/(sqrt(n)*q*G);
	double max_lam;
	max_lam = lam_max * max_tmp / alpha;
	if (lam_min == 0) {
		lambda1.setLinSpaced(nlam - 1, log(max_lam), log(0.0001*max_lam));
		lambda.setLinSpaced(nlam, 0, 0);
		lambda.segment(0, nlam - 1) = lambda1.array().exp().matrix();
	}
	else {
		lambda1.setLinSpaced(nlam, log(max_lam), log(lam_min*max_lam));
		lambda = lambda1.array().exp();
	}
	return lambda;
}
//***-------------------------------------------------------------**
//***-------update the jth row of matrix A with penalty-----------**
MatrixXd updateT4A_penalty(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, 
	MatrixXi beta, VectorXi &activeA, double lambda1)
{
/*
	Input:
	Y is n*q matrix
	Z is n*(K*p*G) matrix
	A is p*r1 matrix
	B is K*r2 matrix
	C is G*r3 matrix
	D is q*r4 matrix
	S is r4*(r1*r2*r3) matrix
	penalty= 1(lasso),2(mcp), and 3(scad)
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	eps is the error to control convergence for outer iteration
	eps1 is the error to control convergence for inner iteration
	max_step is the max step to control convergence for outer iteration
	max_step1 is the max step to control convergence for inner iteration
	nlam is the number of tuning parameters
*/  
	int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, r4 = opts.r4, K = opts.K, G = opts.G, q = opts.q, n = opts.n, p = opts.p;
	int nlam = opts_pen.nlam, gamma = opts_pen.gamma, pen = opts_pen.pen, dfmax = opts_pen.dfmax;
	double alpha = opts_pen.alpha, eps1 = opts.eps1;
	int max_step1 = opts.max_step1, j,jj, active, KG=K*G;

	MatrixXd S1,aj,ajnew,zj;
    VectorXi dims = VectorXi::Constant(4, 0);
    dims[0] = r1; dims[1] = r2; dims[2] = r3; dims[3] = r4;
	S1 = TransferModalUnfoldingsT(S, 4, 1, dims);
	VectorXi id = SEQ(1, p*KG, p);
	MatrixXd Vnew = MatrixXd::Constant(n*q, r1*p, 0);
	MatrixXd Gamma_sqrtn = MatrixXd::Constant(r1, r1*p, 0);
	MatrixXd dcbs = kroneckerProduct(D,kroneckerProduct(C, B))*(S1.transpose());
	MatrixXd A1 = A, V, Gammaj, Gamma_sqrt, V_1, D4,L;
	D4 = D*S*kroneckerProduct(C, kroneckerProduct(B, A)).transpose(); 

	for (j = 0; j < p; j++) {
		V = submatrix_col(Z, id.array() + j)*(dcbs.block(0, 0, KG, r1));
		for (jj = 1; jj < q; jj++) {
			V_1 = submatrix_col(Z, id.array() + j)*(dcbs.block(jj*KG, 0, KG, r1));
			V = rbind_rcpp(V, V_1);
		}
		Gammaj = ((V.transpose()*V).array() / n).matrix();
		L = Gammaj.llt().matrixL();
		Gamma_sqrtn.block(0, j*r1, r1, r1) = UpTriangularInv(L.transpose());
		Vnew.block(0, j*r1, n*q, r1) = QbyR(V, UpTriangularInv(L.transpose()),1);
		A1.row(j) = (QbyR(A.row(j).transpose(), L.transpose(), 0)).transpose();
	}		

	MatrixXd Anew = A1;
	MatrixXd r = Y - Z * (D4.transpose());
	r.resize(n*q, 1);
	VectorXd ajnorm_old, ajnorm;
	ajnorm_old = ajnorm = VectorXd::Constant(p, 0);
	int converged1, step = 0;
	while (step<max_step1) 
	{
		step++;
		active = 0;
		for (j = 0; j < p; j++)
			if (ajnorm[j] != 0) active = active + 1;
		if (active>dfmax) {
			beta = MatrixXi::Constant(p*r1, nlam, -9);
			return A;
		}
		for (j = 0; j < p;j++) {
			aj = Anew.row(j).transpose();
			zj = Vnew.block(0, j*r1, n*q, r1).transpose()*r/n + aj;
			ajnew = updateAj(zj, n, r1, lambda1, alpha, gamma, pen);
			r = r - Vnew.block(0, j*r1, n*q, r1)*(ajnew-aj);
			Anew.row(j) = ajnew.transpose();
			ajnorm[j] = ajnew.norm();
		}
		converged1 = 1;
		for (j = 0; j < p;j++) {
				if (ajnorm[j] != 0 && ajnorm_old[j] != 0) {
					if ((A1.row(j) - Anew.row(j)).norm() / ajnorm_old[j]>eps1) {
						converged1 = 0; break;
					}
				}
				else if (ajnorm[j] == 0 && ajnorm_old[j] != 0) {
					converged1 = 0; break;
				}
				else if (ajnorm[j] != 0 && ajnorm_old[j] == 0) {
					converged1 = 0; break;
				}
			}
		if (converged1) break;
		A1 = Anew;
		ajnorm_old = ajnorm;
    }//end while
	for (j = 0; j<p; j++) {
		Anew.row(j) = (QbyR(Anew.row(j).transpose(),Gamma_sqrtn.block(0, j*r1, r1, r1),0)).transpose();
		if (ajnorm[j]) activeA[j] = 1;
	}
	return Anew;
}

//***-------------------------------------------------------------**
//***------Old main function: Estimation with penalizing functions in a whole column -----------------**
// [[Rcpp::export]]
List EstPenColumnT4(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, 
	VectorXd lambda, List optsList, List optsList_pen)
{
	/*
	Input:
	Y is n*q matrix
	Z is n*(K*p*G) matrix
	A is p*r1 matrix
	B is K*r2 matrix
	C is G*r3 matrix
	D is q*r4 matrix
	S is r4*(r1*r2*r3) matrix
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	gamma is another tuning parameter for MCP and SCAD
	penalty= 1(lasso),2(mcp), and 3(scad)
	dfmax is the preset maximum digrees of freedom
	eps is the error to control convergence for outer iteration
	eps1 is the error to control convergence for inner iteration
	max_step is the max step to control convergence for outer iteration
	max_step1 is the max step to control convergence for inner iteration
	nlam is the number of tuning parameters

	Output:
	Dnew is a estimator of D(4) = D%*%S%*%kronecker(t(C),kronecker(t(B),t(A)))
	*/
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<double>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);
	opts.max_step1 = as<double>(optsList["max_step1"]);	
	opts.r1 = as<int>(optsList["r1"]);
	opts.r2 = as<int>(optsList["r2"]);
	opts.r3 = as<int>(optsList["r3"]);
	opts.r4 = as<int>(optsList["r4"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.G = as<int>(optsList["G"]);
	opts.degr = as<int>(optsList["degr"]);
	opts.K = as<int>(optsList["K"]);
//    opts.n = Y.cols();
    opts.n = Y.rows();
		
	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma = as<double>(optsList_pen["gamma"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.size();	
	
	
	int l,j, step, max_step=opts.max_step, nlam=opts_pen.nlam, G = opts.G, K = opts.K, KG = K*G, p = opts.p;
	int q = opts.q, n = opts.n, r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, r4 = opts.r4;
	double  likhd0 = pow(10, 6), lambda1, likhd1, eps = opts.eps;
	MatrixXd Dnew, Anew, Bnew, Snew, Cnew, Dn, Z1, A1, Z2, A2;
	VectorXi activeA = VectorXi::Constant(p, 1), convergence1 = VectorXi::Constant(5, 1);	
	VectorXi df = VectorXi::Constant(nlam, 0);
	MatrixXi betapath = MatrixXi::Constant(p, nlam, 0);
	VectorXd likhd = VectorXd::Constant(nlam, 0);
	MatrixXd Apath, Bpath, Cpath, Dpath, Spath, temp;
	Apath = MatrixXd::Constant(p*r1, nlam, 0);
	Bpath = MatrixXd::Constant(K*r2, nlam, 0);
	Cpath = MatrixXd::Constant(G*r3, nlam, 0);
	Dpath = MatrixXd::Constant(q*r4, nlam, 0);
	Spath = MatrixXd::Constant(r1*r1*r3*r4, nlam, 0);	
    Anew = A;
    A1 = A;
	Z1 = Z;
//    Rcout << "nlam= " << nlam << std::endl;
	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step = 0;
		while (step<max_step) {
			step ++;
			Snew = updateT4S(Y, Z1, A1, B, C, D);
			Dn = D * Snew * kroneckerProduct(C,kroneckerProduct(B, A1)).transpose();
			likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            Rcout << "step= " << step << "        S likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
			if (likhd1<likhd0) {
				S = Snew;
				likhd0 = likhd1;
			}
			else convergence1[0]=0;

			Dnew = updateT4D(Y, Z1, S, A1, B, C, D);
			Dn = Dnew * S*kroneckerProduct(C,kroneckerProduct(B, A1)).transpose();
			likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();	
            Rcout << "step= " << step << "        D likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
			if (likhd1<likhd0) {
				D = Dnew;
				likhd0 = likhd1;
			}
			else convergence1[1]=0;	

			Cnew = updateT4C(Y, Z1, S, A1, B, C, D);
			Dn = D * S * kroneckerProduct(Cnew,kroneckerProduct(B, A1)).transpose();
			likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            Rcout << "step= " << step << "        C likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
			if (likhd1<likhd0) {
				C = Cnew;
				likhd0 = likhd1;
			}
			else convergence1[2]=0;

			Bnew = updateT4B(Y, Z1, S, A1, B, C, D);
			Dn = D * S * kroneckerProduct(C,kroneckerProduct(Bnew, A1)).transpose();
			likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            Rcout << "step= " << step << "        B likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
			//stop("stop here!");
			if (likhd1<likhd0) {
				B = Bnew;
				likhd0 = likhd1;	
			}
			else convergence1[3]=0;

			for(j=0;j<p;j++) activeA[j] = 0;
			Anew = updateT4A_penalty(Y, Z, S, A, B, C, D, betapath, activeA, lambda1);
//            Rcout << "activeA.sum()= " << activeA.sum() << " r1= " << r1 << std::endl;
			if(activeA.sum()<r1){
				A1 = A;
				Z1 = Z;
				break;
			}
			else{
				Z2 = extractColsZ(Z,p,KG,activeA);
				A2 = extractRows(Anew, activeA);
				Dn = D * S * kroneckerProduct(C,kroneckerProduct(B, A2)).transpose();
				likhd1 = (Y - Z2 * Dn.transpose()).squaredNorm();
                Rcout << "likhd1= " << likhd1 << " likhd0= " << likhd0 << std::endl;
				if (likhd1<likhd0) {
					A = Anew;
					Z1 = Z2;
					A1 = A2;
					if ((likhd0 - likhd1) / likhd0<eps) break;
					else  likhd0 = likhd1;	
				}
				else convergence1[4]=0;
			}
            Rcout << "step= " << step << "        A likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
            
//            if (step==5) stop("stop for test!");
			if(convergence1.sum()==0) break;
			else for(j=0;j<5;j++) convergence1[j] = 1;
		} //end while		
		for(j=0;j<p;j++) activeA[j] = 0;
		for(j=0;j<p;j++) if(A.row(j).norm()) activeA[j] = 1;
		df[l] = activeA.sum();
		likhd[l] = likhd0;
		betapath.col(l) = activeA;

		temp = A; temp.resize(p*r1, 1);		Apath.col(l) = temp;
		temp = B; temp.resize(K*r2, 1);		Bpath.col(l) = temp;
		temp = C; temp.resize(G*r3, 1);		Cpath.col(l) = temp;
		temp = D; temp.resize(q*r4, 1);		Dpath.col(l) = temp;
		temp = S; temp.resize(r1*r1*r3*r4, 1);		Spath.col(l) = temp;
	}// end for
	return List::create(Named("likhd") = likhd, Named("betapath") = betapath, Named("df") = df, Named("lambda")=lambda, 
	                    Named("Spath") = Spath, Named("Apath") = Apath, Named("Bpath") = Bpath, Named("Cpath") = Cpath, Named("Dpath") = Dpath);
}

//***--------------------------------------------------------------**
//***----------- Estimation with penalizing functions in a whole column by CV---------------------**
// [[Rcpp::export]]
List EstPenColumnT4CV(MatrixXd Y, MatrixXd Z, MatrixXd Ytest, MatrixXd Ztest, MatrixXd S, MatrixXd A, MatrixXd B, 
                    MatrixXd C, MatrixXd D, VectorXd lambda, List optsList, List optsList_pen)
{
	/*
	Input:
	Y is n*q matrix
	Z is n*(K*p*G) matrix
	A is p*r1 matrix
	B is K*r2 matrix
	C is G*r3 matrix
	D is q*r4 matrix
	S is r4*(r1*r2*r3) matrix
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	gamma is another tuning parameter for MCP and SCAD
	penalty= 1(lasso),2(mcp), and 3(scad)
	dfmax is the preset maximum digrees of freedom
	eps is the error to control convergence for outer iteration
	eps1 is the error to control convergence for inner iteration
	max_step is the max step to control convergence for outer iteration
	max_step1 is the max step to control convergence for inner iteration
	nlam is the number of tuning parameters

	Output:
	Dnew is a estimator of D(4) = D%*%S%*%kronecker(t(C),kronecker(t(B),t(A)))
	*/
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<double>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);
	opts.max_step1 = as<double>(optsList["max_step1"]);	
	opts.r1 = as<int>(optsList["r1"]);
	opts.r2 = as<int>(optsList["r2"]);
	opts.r3 = as<int>(optsList["r3"]);
	opts.r4 = as<int>(optsList["r4"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.G = as<int>(optsList["G"]);
	opts.degr = as<int>(optsList["degr"]);
	opts.K = as<int>(optsList["K"]);
//    opts.n = Y.cols();
    opts.n = Y.rows();
		
	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma = as<double>(optsList_pen["gamma"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.size();	
	
	int l,j, step, max_step=opts.max_step, nlam=opts_pen.nlam, KG = opts.K*opts.G, p = opts.p, r1 = opts.r1, n = opts.n;
	double  likhd0 = pow(10, 6), lambda1, likhd1, eps = opts.eps;	
	MatrixXd Dnew, Anew, Bnew, Snew, Cnew, Dn, Z1, A1, Z2, A2;
	VectorXi activeA = VectorXi::Constant(p, 0), convergence1 = VectorXi::Constant(5, 1);
	VectorXi df = VectorXi::Constant(nlam, 0);
	MatrixXi betapath = MatrixXi::Constant(p, nlam, 0);
	VectorXd likhd = VectorXd::Constant(nlam, 0);
    Anew = A;
    A1 = A;
	Z1 = Z;
	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step = 0;
		while (step<max_step) {
			step ++;
			Snew = updateT4S(Y, Z1, A1, B, C, D);
			Dn = D * Snew * kroneckerProduct(C,kroneckerProduct(B, A1)).transpose();
			likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
			//Rcout << "step= " << step << "        S likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;			
			if (likhd1<likhd0) {
				S = Snew;
				likhd0 = likhd1;
			}
			else convergence1[0]=0;
			
			Dnew = updateT4D(Y, Z1, S, A1, B, C, D);
			Dn = Dnew * S*kroneckerProduct(C,kroneckerProduct(B, A1)).transpose();
			likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();	
			//Rcout << "step= " << step << "        D likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
			if (likhd1<likhd0) {
				D = Dnew;
				likhd0 = likhd1;
			}
			else convergence1[1]=0;	

			Cnew = updateT4C(Y, Z1, S, A1, B, C, D);
			Dn = D * S * kroneckerProduct(Cnew,kroneckerProduct(B, A1)).transpose();
			likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
			//Rcout << "step= " << step << "        C likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;			
			if (likhd1<likhd0) {
				C = Cnew;
				likhd0 = likhd1;
			}
			else convergence1[2]=0;

			Bnew = updateT4B(Y, Z1, S, A1, B, C, D);
			Dn = D * S * kroneckerProduct(C,kroneckerProduct(Bnew, A1)).transpose();
			likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
			//Rcout << "step= " << step << "        B likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
			//if(step==3) stop("stop here!");			
			if (likhd1<likhd0) {
				B = Bnew;
				likhd0 = likhd1;	
			}
			else convergence1[3]=0;

			for(j=0;j<p;j++) activeA[j] = 0;
			Anew = updateT4A_penalty(Y, Z, S, A, B, C, D, betapath, activeA, lambda1);
			if(activeA.sum()<r1){
				A1 = A;
				Z1 = Z;
				break;
			}
			else{
				Z2 = extractColsZ(Z,p,KG,activeA);
				A2 = extractRows(Anew, activeA);
				Dn = D * S * kroneckerProduct(C,kroneckerProduct(B, A2)).transpose();
				likhd1 = (Y - Z2 * Dn.transpose()).squaredNorm();
			//Rcout << "step= " << step << "        A likhd1= " << likhd1<< "        likhd0= " << likhd0<< std::endl;
			//Rcout << "A likhd1= " << (likhd0 - likhd1) / likhd0<< "        eps= " << eps<< std::endl;
			//stop("stop here!");				
				if (likhd1<likhd0) {
					A = Anew;
					Z1 = Z2;
					A1 = A2;
					if ((likhd0 - likhd1) / likhd0<eps) break;
					else  likhd0 = likhd1;						
				}
				else convergence1[4]=0;
			}		
			//Rcout << "convergence1.sum()= " << convergence1.sum() << std::endl;
			if(convergence1.sum()==0) break;
			else for(j=0;j<5;j++) convergence1[j] = 1;
			//Rcout << "Z1.cols= " << Z1.cols()<< "        A.rows= " << A.rows()<< "  activeA.sum = "<< activeA.sum()<< std::endl;
			//if(step==2) stop("stop here!");
		} //end while		
		//stop("stop here!");
		for(j=0;j<p;j++) activeA[j] = 0;
		for(j=0;j<p;j++) if(A.row(j).norm()) activeA[j] = 1;
		df[l] = activeA.sum();
		Z2 = extractColsZ(Ztest,p,KG,activeA);
		A2 = extractRows(A, activeA);
		Dn = D * S * kroneckerProduct(C,kroneckerProduct(B, A2)).transpose();
		likhd[l] = (Ytest - Z2 * Dnew.transpose()).squaredNorm();
		betapath.col(l) = activeA;		
	}// end for
	return List::create(Named("likhd") = likhd, Named("df") = df, Named("betapath")=betapath);
}

