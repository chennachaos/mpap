#ifndef incl_EIGEN_COND_NUMBER_H
#define incl_EIGEN_COND_NUMBER_H


#include <Eigen/Core>
#include <Eigen/Sparse>

#include <cmath>
#include <random>
#include <vector>

using std::cout;
static std::mt19937 engine;
static const int itmax=10;

namespace myMatrix
{

class OneNormEst
{
public:
    OneNormEst(const int& n, const int& t=4): m_n(n), m_t(t) {}
    ~OneNormEst() {}

    template<typename A_operator, typename A_transpose_operator>
    bool ANorm(const A_operator& applyA, A_transpose_operator& applyA_trans, double& norm)
    {
        Eigen::MatrixXd X(m_n, m_t), Y(m_n, m_t), Z(m_n, m_t);

        std::uniform_int_distribution<int> rand_gen(0, m_n-1);
        auto random_plus_minus_1_func = [&](const double& ) -> double
        {
            if (rand_gen(engine)>m_n/2)
                return 1;
            else
                return -1;
        };

        X = X.unaryExpr(random_plus_minus_1_func);
        X.col(0).setOnes();
        X /= m_n;

        Eigen::MatrixXd S(m_n, m_t), S_old(m_n, m_t); // sign matrix
        Eigen::MatrixXi prodS(m_t, m_t); 
        S.setZero(); 
        
        double est=0., est_old = 0.;

        int ind_best;
        std::vector<int> indices(m_n);
        std::vector<int> ind_hist;

        for (int k(0); k<itmax; ++k)
        {
            //applyA(X, Y); // Y = A * X
            Y = applyA * X;
        
            int ind(-1); 
            for(int i(0); i<m_t; ++i)
            {
                double norm = Y.col(i).cwiseAbs().sum(); // norm = {||Y||}_1 
                if (norm > est) {
                    est = norm;
                    ind = indices[i];
                }
            }

            if (est > est_old || k==1)
                ind_best = ind;
        
            if (est < est_old && k >= 1)
            {
                norm = est_old;
                return true; 
            }
        
            est_old = est;
            S_old = S;
        
            // S = sign(Y)
            S = Y.unaryExpr([&] (const double& coeff) -> double
            {
                if(coeff >= 0.) return 1.;
                else            return -1.;
            });

            prodS = (S_old.transpose() * S).matrix().cast<int>();

            // if all cols are parallel, prodS will have all entries = n
            if (prodS.cwiseAbs().colwise().sum().sum() == m_n*m_n*m_t) {
                norm = est;
                return true;
            }
        
            // if S_old(i) is parallel to S(j), replace S(j) with random 
            for(int i(0); i<m_t; ++i)
            {
                for(int j(0); j<m_t; ++j)
                {
                    if(prodS.coeff(i, j)==m_n) 
                        S.col(j) = S.col(j).unaryExpr(random_plus_minus_1_func);
                }
            }

            // if S(i) is parallel to S(j), replace S(j) with random 
            prodS = (S.transpose()*S).matrix().cast<int>();

            for(int i(0); i<m_t; ++i)
            {
                for(int j(i+1); j<m_t; ++j)
                {
                    if(prodS.coeff(i, j)==m_n) 
                        S.col(j) = S.col(j).unaryExpr(random_plus_minus_1_func);
                }
            }

            //applyA_trans(S, Z);
            Z = applyA_trans * S;

            Eigen::VectorXd h = Z.cwiseAbs().rowwise().maxCoeff();
            
            if(k >= 1 && h.maxCoeff() == h.coeff(ind_best) )
            {
              norm = est;
              return true;
            }
        
            for(int i(0); i<m_n; ++i)
              indices[i] = i;
        
            // sort indices by entries in h
            //std::sort(indices.begin(), indices.end(), [&](int& left, int& right)
            //{ return h.coeff(left) > h.coeff(right);  });

            std::sort(indices.begin(), indices.end(), [&](const int left, const int right)
            { return h.coeff(left) > h.coeff(right);  });

            int n_found(0);
            for(auto i=indices.begin(); i!=indices.begin()+m_t; ++i) // is indices contained in ind_hist?
                if(std::find(ind_hist.begin(), ind_hist.end(), *i) != ind_hist.end())
                    ++n_found;

            if(n_found == m_t) {
                norm = est;
                return true;
            }

            std::vector<int> new_indices;

            if(k>0) 
            {
                new_indices.reserve(m_t);
                int count(0);
                for(auto it=indices.begin()+m_t; it !=indices.end() && count <m_t; ++it)
                {
                    if(std::find(ind_hist.begin(), ind_hist.end(), *it) == ind_hist.end()) {
                        new_indices.push_back(*it); ++count;
                    }
                }
                new_indices.swap(indices);
            }
            assert(indices.size()>0);

            X.setZero();
            for(auto i=0; i<m_t; ++i)
            {    
                X.coeffRef(indices[i], i) = 1; //Set X(:, j) = e_ind-j 
            
                ind_hist.push_back(indices[i]);
            }
        }
        norm = est;

        return true;

    }

private:
    int m_n; // rows
    int m_t; // cols

};

}


#endif