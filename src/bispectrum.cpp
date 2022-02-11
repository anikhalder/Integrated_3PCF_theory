#include <bispectrum.hpp>
#include <cosmology_utils.h>
#include <integration_utils.h>
#include <constants.h>
#include <assert.h>

namespace
{
    // Scoccimarro and Couchman bispectrum fitting formula parameter values

    const double a1_SC =  0.25;
    const double a2_SC =  3.5;
    const double a3_SC =  2.0;
    const double a4_SC =  1.0;
    const double a5_SC =  2.0;
    const double a6_SC = -0.2;

    // Gil-Marin bispectrum fitting formula parameter values

    const double a1_GM =  0.484;
    const double a2_GM =  3.740;
    const double a3_GM = -0.849;
    const double a4_GM =  0.392;
    const double a5_GM =  1.013;
    const double a6_GM = -0.575;
    const double a7_GM =  0.128;
    const double a8_GM = -0.722;
    const double a9_GM = -0.926;

    const double f_sq = 7.0; // squeezed factor for response function approach
}

bool is_triangle_closed(const double &k_1, const double &k_2, const double &k_3)
{
    double k[4], kt;
    k[1]=k_1, k[2]=k_2, k[3]=k_3;
    for(int i=1;i<=3;i++) // sorting k[i] such that k[1]>=k[2]>=k[3]
    {
        for(int j=i+1;j<=3;j++)
        {
            if(k[i]<k[j])
            {
                kt=k[j];
                k[j]=k[i];
                k[i]=kt;
            }
        }
    }

    if(k[1]>k[2]+k[3])
    {
        //printf("Error: triangle is not closed! \n");
        //printf("k1, k2, k3 = %f, %f, %f \n", k[1], k[2], k[3]);
        return false;
    }
    else
        return true;
}


// ######################################################################################

// bispectrum prescription

double B(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    if (z > class_obj->get_z_max_pk())
        return 0.0;

    double k_max_pk = class_obj->get_k_max_pk();
    if (k_1 > k_max_pk || k_2 > k_max_pk || k_3 > k_max_pk)
        return 0.0;

    if (k_1 == 0 || k_2 == 0 || k_3 == 0)
        return 0.0;

    if (apply_T17_corrections && (T17_box_correction(k_1,z,class_obj) == 0 || T17_box_correction(k_2,z,class_obj) == 0 || T17_box_correction(k_3,z,class_obj) == 0))
        return 0.0;

    double k[4], kt;
    k[1]=k_1, k[2]=k_2, k[3]=k_3;
    for(int i=1;i<=3;i++) // sorting k[i] such that k[1]>=k[2]>=k[3]
    {
        for(int j=i+1;j<=3;j++)
        {
            if(k[i]<k[j])
            {
                kt=k[j];
                k[j]=k[i];
                k[i]=kt;
            }
        }
    }

    if(k[1]>k[2]+k[3])
    {
        //printf("Error: triangle is not closed! \n");
        return 0.0;
    }

    double k_h = k[1]; // hard mode --> MAX(k1,k2,k3)
    double k_m = k[2];
    double k_s = k[3]; // soft mode --> MIN(k1,k2,k3)

    //double B = B_tree(k_h, k_m, k_s, z, class_obj, false);
    //double B = B_1_loop_hcubature(k_h, k_m, k_s, z, class_obj, false);
    //double B = B_SC(k_h, k_m, k_s, z, class_obj, true);
    //double B = B_GM(k_h, k_m, k_s, z, class_obj, true);
    //double B = B_bihalofit(k_h, k_m, k_s, z, class_obj, false);
    //double B = B_squeezed_RF(k_h, k_m, k_s, z, class_obj, true);

    // Use hybrid model of bispectra for squeezed and non-squeezed triangle configurations

    
    double B;

    if (k_m < f_sq*k_s)
    {
        // non-squeezed configurations

        //B = 0.0;
        //B = B_tree(k_h, k_m, k_s, z, class_obj, false);
        //B = B_1_loop_hcubature(k_h, k_m, k_s, z, class_obj, false);
        //B = B_SC(k_h, k_m, k_s, z, class_obj, true);
        B = B_GM(k_h, k_m, k_s, z, class_obj, true);
        //B = B_bihalofit(k_h, k_m, k_s, z, class_obj, false);
        //B = B_squeezed_RF(k_h, k_m, k_s, z, class_obj, true);
    }

    else
    {
        // squeezed configurations

        //B = 0.0;
        //B = B_tree(k_h, k_m, k_s, z, class_obj, false);
        //B = B_1_loop_hcubature(k_h, k_m, k_s, z, class_obj, false);
        //B = B_SC(k_h, k_m, k_s, z, class_obj, true);
        //B = B_GM(k_h, k_m, k_s, z, class_obj, true);
        //B = B_bihalofit(k_h, k_m, k_s, z, class_obj, false);
        B = B_squeezed_RF(k_h, k_m, k_s, z, class_obj, true);

    }

    if (apply_T17_corrections)
        B *= T17_shell_correction(k_h,class_obj)*T17_shell_correction(k_m,class_obj)*T17_shell_correction(k_s,class_obj);

    return B;
}

// ######################################################################################

// tree-level SPT matter bispectrum

double B_tree(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    // Already ensure beforehand that the input k_i modes form a closed triangle

//    if (class_obj->get_f_NL_local() != 0)
//        B += B_primordial_local(k_1,k_2,k_3,z,class_obj);

//    else if (class_obj->get_f_NL_equilateral() != 0)
//        B += B_primordial_equilateral(k_1,k_2,k_3,z,class_obj);

//    else if (class_obj->get_f_NL_orthogonal() != 0)
//        B += B_primordial_orthogonal(k_1,k_2,k_3,z,class_obj);

    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
    double Pk_1 = class_obj->pk(k_1,z,use_pk_nl);
    double Pk_2 = class_obj->pk(k_2,z,use_pk_nl);
    double Pk_3 = class_obj->pk(k_3,z,use_pk_nl);

    double B = 2.0*(F2_EdS(k_1,k_2,k_3)*Pk_1*Pk_2
                   +F2_EdS(k_2,k_3,k_1)*Pk_2*Pk_3
                   +F2_EdS(k_3,k_1,k_2)*Pk_3*Pk_1);

    return B;
}


// ######################################################################################

// components for tracer bispectrum

double B_P(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    if (z > class_obj->get_z_max_pk())
        return 0.0;

    double k_max_pk = class_obj->get_k_max_pk();
    if (k_1 > k_max_pk || k_2 > k_max_pk || k_3 > k_max_pk)
        return 0.0;

    if (k_1 == 0 || k_2 == 0 || k_3 == 0)
        return 0.0;

    if (apply_T17_corrections && (T17_box_correction(k_1,z,class_obj) == 0 || T17_box_correction(k_2,z,class_obj) == 0 || T17_box_correction(k_3,z,class_obj) == 0))
        return 0.0;

    double k[4], kt;
    k[1]=k_1, k[2]=k_2, k[3]=k_3;
    for(int i=1;i<=3;i++) // sorting k[i] such that k[1]>=k[2]>=k[3]
    {
        for(int j=i+1;j<=3;j++)
        {
            if(k[i]<k[j])
            {
                kt=k[j];
                k[j]=k[i];
                k[i]=kt;
            }
        }
    }

    if(k[1]>k[2]+k[3])
    {
        //printf("Error: triangle is not closed! \n");
        return 0.0;
    }

    double k_h = k[1]; // hard mode --> MAX(k1,k2,k3)
    double k_m = k[2];
    double k_s = k[3]; // soft mode --> MIN(k1,k2,k3)

    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
    double Pk_h = class_obj->pk(k_h,z,use_pk_nl);
    double Pk_m = class_obj->pk(k_m,z,use_pk_nl);
    double Pk_s = class_obj->pk(k_s,z,use_pk_nl);

    double B_P = Pk_h + Pk_m + Pk_s;

    return B_P;
}

double B_PP(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    if (z > class_obj->get_z_max_pk())
        return 0.0;

    double k_max_pk = class_obj->get_k_max_pk();
    if (k_1 > k_max_pk || k_2 > k_max_pk || k_3 > k_max_pk)
        return 0.0;

    if (k_1 == 0 || k_2 == 0 || k_3 == 0)
        return 0.0;

    if (apply_T17_corrections && (T17_box_correction(k_1,z,class_obj) == 0 || T17_box_correction(k_2,z,class_obj) == 0 || T17_box_correction(k_3,z,class_obj) == 0))
        return 0.0;

    double k[4], kt;
    k[1]=k_1, k[2]=k_2, k[3]=k_3;
    for(int i=1;i<=3;i++) // sorting k[i] such that k[1]>=k[2]>=k[3]
    {
        for(int j=i+1;j<=3;j++)
        {
            if(k[i]<k[j])
            {
                kt=k[j];
                k[j]=k[i];
                k[i]=kt;
            }
        }
    }

    if(k[1]>k[2]+k[3])
    {
        //printf("Error: triangle is not closed! \n");
        return 0.0;
    }

    double k_h = k[1]; // hard mode --> MAX(k1,k2,k3)
    double k_m = k[2];
    double k_s = k[3]; // soft mode --> MIN(k1,k2,k3)

    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
    double Pk_h = class_obj->pk(k_h,z,use_pk_nl);
    double Pk_m = class_obj->pk(k_m,z,use_pk_nl);
    double Pk_s = class_obj->pk(k_s,z,use_pk_nl);

    double B_PP = Pk_h*Pk_m + Pk_h*Pk_s + Pk_m*Pk_s;

    return B_PP;
}

double B_S2PP(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    if (z > class_obj->get_z_max_pk())
        return 0.0;

    double k_max_pk = class_obj->get_k_max_pk();
    if (k_1 > k_max_pk || k_2 > k_max_pk || k_3 > k_max_pk)
        return 0.0;

    if (k_1 == 0 || k_2 == 0 || k_3 == 0)
        return 0.0;

    if (apply_T17_corrections && (T17_box_correction(k_1,z,class_obj) == 0 || T17_box_correction(k_2,z,class_obj) == 0 || T17_box_correction(k_3,z,class_obj) == 0))
        return 0.0;

    double k[4], kt;
    k[1]=k_1, k[2]=k_2, k[3]=k_3;
    for(int i=1;i<=3;i++) // sorting k[i] such that k[1]>=k[2]>=k[3]
    {
        for(int j=i+1;j<=3;j++)
        {
            if(k[i]<k[j])
            {
                kt=k[j];
                k[j]=k[i];
                k[i]=kt;
            }
        }
    }

    if(k[1]>k[2]+k[3])
    {
        //printf("Error: triangle is not closed! \n");
        return 0.0;
    }

    double k_h = k[1]; // hard mode --> MAX(k1,k2,k3)
    double k_m = k[2];
    double k_s = k[3]; // soft mode --> MIN(k1,k2,k3)

    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
    double Pk_h = class_obj->pk(k_h,z,use_pk_nl);
    double Pk_m = class_obj->pk(k_m,z,use_pk_nl);
    double Pk_s = class_obj->pk(k_s,z,use_pk_nl);

    double S2_h_m = S2(k_h,k_m,k_s);
    double S2_h_s = S2(k_h,k_s,k_m);
    double S2_m_s = S2(k_m,k_s,k_h);

    double B_S2PP = S2_h_m*Pk_h*Pk_m + S2_h_s*Pk_h*Pk_s + S2_m_s*Pk_m*Pk_s;

    return B_S2PP;
}

double B_PaPb(const double &k_a, const double &k_b, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    if (z > class_obj->get_z_max_pk())
        return 0.0;

    double k_max_pk = class_obj->get_k_max_pk();
    if (k_a > k_max_pk || k_b > k_max_pk)
        return 0.0;

    if (k_a == 0 || k_b == 0)
        return 0.0;

    if (apply_T17_corrections && (T17_box_correction(k_a,z,class_obj) == 0 || T17_box_correction(k_b,z,class_obj) == 0))
        return 0.0;

    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
    double Pk_a = class_obj->pk(k_a,z,use_pk_nl);
    double Pk_b = class_obj->pk(k_b,z,use_pk_nl);

    double B_PaPb = Pk_a*Pk_b;

    return B_PaPb;
}

double B_S2PaPb(const double &k_a, const double &k_b, const double &k_c, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    if (z > class_obj->get_z_max_pk())
        return 0.0;

    double k_max_pk = class_obj->get_k_max_pk();
    if (k_a > k_max_pk || k_b > k_max_pk || k_c > k_max_pk)
        return 0.0;

    if (k_a == 0 || k_b == 0 || k_c == 0)
        return 0.0;

    if (apply_T17_corrections && (T17_box_correction(k_a,z,class_obj) == 0 || T17_box_correction(k_b,z,class_obj) == 0 || T17_box_correction(k_c,z,class_obj) == 0))
        return 0.0;

    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
    double Pk_a = class_obj->pk(k_a,z,use_pk_nl);
    double Pk_b = class_obj->pk(k_b,z,use_pk_nl);

    double S2_a_b = S2(k_a,k_b,k_c);

    double B_S2PaPb = S2_a_b*Pk_a*Pk_b;

    return B_S2PaPb;
}


// ######################################################################################

// 1-loop SPT matter bispectrum

double evaluate_B_1_loop_integrand(const std::vector<double> &q_1_vec, const std::vector<double> &k_a_vec, const std::vector<double> &k_b_vec,
                                   const std::vector<double> &k_c_vec, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    if (!is_vec_alright(q_1_vec))
        std::cout << "NaN encountered in q_1_vec!" << std::endl;

    if (!is_vec_alright(k_a_vec))
        std::cout << "NaN encountered in k_a_vec!" << std::endl;

    if (!is_vec_alright(k_b_vec))
        std::cout << "NaN encountered in k_b_vec!" << std::endl;

    if (!is_vec_alright(k_c_vec))
        std::cout << "NaN encountered in k_c_vec!" << std::endl;

    std::vector<double> pk_A_vec = k_a_vec; // p or plus
    std::vector<double> pk_B_vec = k_b_vec;
    std::vector<double> pk_C_vec = k_c_vec;
    std::vector<double> pq_1_vec = q_1_vec;

    std::vector<double> mk_A_vec = q_mA_vec(k_a_vec); // m for minus
    std::vector<double> mk_B_vec = q_mA_vec(k_b_vec);
    std::vector<double> mk_C_vec = q_mA_vec(k_c_vec);
    std::vector<double> mq_1_vec = q_mA_vec(q_1_vec);

    std::vector<double> pAm1_vec = q_ApB_vec(pk_A_vec,mq_1_vec);
    std::vector<double> mAp1_vec = q_ApB_vec(mk_A_vec,pq_1_vec);

    std::vector<double> pBp1_vec = q_ApB_vec(pk_B_vec,pq_1_vec);
    std::vector<double> mBm1_vec = q_ApB_vec(mk_B_vec,mq_1_vec);

    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
    double P_A = class_obj->pk(k_a_vec.at(3),z,false);
    double P_B = class_obj->pk(k_b_vec.at(3),z,false);
    double P_C = class_obj->pk(k_c_vec.at(3),z,false);

    double P_1 = class_obj->pk(q_1_vec.at(3),z,false);
    double P_pAm1 = class_obj->pk(pAm1_vec.at(3),z,false);
    double P_pBp1 = class_obj->pk(pBp1_vec.at(3),z,false);

    double B_222_integrand = 8*P_1*F2_sym(pq_1_vec,pAm1_vec)*F2_sym(mq_1_vec,pBp1_vec)*F2_sym(mBm1_vec,mAp1_vec)*P_pAm1*P_pBp1;

    if (isnan(B_222_integrand))
    {
        std::cout << "NaN encountered in B_222_integrand!" << std::endl;
        std::cout << k_b_vec.at(3) << " " << q_1_vec.at(3) << " " << pBp1_vec.at(3) << std::endl;
        std::cout << P_1 << " " << F2_sym(pq_1_vec,pAm1_vec) << " " << F2_sym(mq_1_vec,pBp1_vec) << " " <<  F2_sym(mBm1_vec,mAp1_vec) << " " <<
                     P_pAm1 << " " <<  P_pBp1 << std::endl;
    }

    std::vector<double> pBm1_vec = q_ApB_vec(pk_B_vec,mq_1_vec);
    std::vector<double> mBp1_vec = q_ApB_vec(mk_B_vec,pq_1_vec);

    std::vector<double> pCm1_vec = q_ApB_vec(pk_C_vec,mq_1_vec);
    std::vector<double> mCp1_vec = q_ApB_vec(mk_C_vec,pq_1_vec);

    double P_pBm1 = class_obj->pk(pBm1_vec.at(3),z,false);
    double P_pCm1 = class_obj->pk(pCm1_vec.at(3),z,false);

    double B_123_I_integrand = 6*P_1*(P_A*F2_sym(pq_1_vec,pBm1_vec)*F3_sym(mk_A_vec,mq_1_vec,mBp1_vec)*P_pBm1+
                                      P_A*F2_sym(pq_1_vec,pCm1_vec)*F3_sym(mk_A_vec,mq_1_vec,mCp1_vec)*P_pCm1+
                                      P_B*F2_sym(pq_1_vec,pAm1_vec)*F3_sym(mk_B_vec,mq_1_vec,mAp1_vec)*P_pAm1+
                                      P_B*F2_sym(pq_1_vec,pCm1_vec)*F3_sym(mk_B_vec,mq_1_vec,mCp1_vec)*P_pCm1+
                                      P_C*F2_sym(pq_1_vec,pAm1_vec)*F3_sym(mk_C_vec,mq_1_vec,mAp1_vec)*P_pAm1+
                                      P_C*F2_sym(pq_1_vec,pBm1_vec)*F3_sym(mk_C_vec,mq_1_vec,mBp1_vec)*P_pBm1);

    if (isnan(B_123_I_integrand))
    {
        std::cout << "NaN encountered in B_123_I_integrand!" << std::endl;

    }

    double F2_A_B = F2_sym(pk_A_vec,pk_B_vec);
    double F2_A_C = F2_sym(pk_A_vec,pk_C_vec);
    double F2_B_C = F2_sym(pk_B_vec,pk_C_vec);

    mq_1_vec = q_mA_perturbed_vec (q_1_vec); // perturb slightly to avoid 0 in denominator
    double F3_A_1_m1 = F3_sym(pk_A_vec,pq_1_vec,mq_1_vec);
    double F3_B_1_m1 = F3_sym(pk_B_vec,pq_1_vec,mq_1_vec);
    double F3_C_1_m1 = F3_sym(pk_C_vec,pq_1_vec,mq_1_vec);

    double B_123_II_integrand = 6*P_1*(F2_A_B*P_A*P_B*F3_B_1_m1+
                                       F2_A_C*P_A*P_C*F3_C_1_m1+
                                       F2_A_B*P_A*P_B*F3_A_1_m1+
                                       F2_B_C*P_B*P_C*F3_C_1_m1+
                                       F2_A_C*P_A*P_C*F3_A_1_m1+
                                       F2_B_C*P_B*P_C*F3_B_1_m1);

    if (isnan(B_123_II_integrand))
        std::cout << "NaN encountered in B_123_II_integrand!" << std::endl;

    double B_114_integrand = 12*P_1*(P_A*P_B*F4_sym(mk_A_vec,mk_B_vec,mq_1_vec,pq_1_vec)+
                                     P_A*P_C*F4_sym(mk_A_vec,mk_C_vec,mq_1_vec,pq_1_vec)+
                                     P_B*P_C*F4_sym(mk_B_vec,mk_C_vec,mq_1_vec,pq_1_vec));

    if (isnan(B_114_integrand))
        std::cout << "NaN encountered in B_114_integrand!" << std::endl;

    return B_222_integrand + B_123_I_integrand + B_123_II_integrand + B_114_integrand;
}

// h-cubature

int B_1_loop_hcubature_integrand(unsigned ndim, const double *q, void *params, unsigned fdim, double *value)
{
    assert(ndim == 3);
    assert(fdim == 1);

    struct params_B_1_loop_integrand *p = static_cast<params_B_1_loop_integrand *>(params);

    double r = q[0];
    double cos_theta = q[1];
    double phi = q[2];
    double sin_theta = sqrt(1-cos_theta*cos_theta);

    std::vector<double> q_1_vec = {r*sin_theta*cos(phi),r*sin_theta*sin(phi),r*cos_theta,r};

    value[0] = r*r*evaluate_B_1_loop_integrand(q_1_vec, p->k_1_vec, p->k_2_vec, p->k_3_vec, p->z, p->class_obj, p->use_pk_nl);

    return 0;
}

double B_1_loop_hcubature(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    // Already ensure beforehand that the input k_i modes form a closed triangle

    double k_max_integral = 2.0*class_obj->get_h();
    if (k_1 > k_max_integral && k_2 > k_max_integral && k_3 > k_max_integral)
        return 0.0;

    double k_min_integral = 0.001*class_obj->get_h();
    if (k_1 < k_min_integral && k_2 < k_min_integral && k_3 < k_min_integral)
        return B_tree(k_1, k_2, k_3, z, class_obj, use_pk_nl);

    double result = 0;              // the result from the integration
    double error = 0;               // the estimated error from the integration

    double cos_phi_12 = 0.5*(k_3*k_3-k_1*k_1-k_2*k_2)/(k_1*k_2); // law of cosines
    double sin_phi_12 = sqrt(1 - cos_phi_12*cos_phi_12);

    if (isnan(sin_phi_12))
        sin_phi_12 = 0;

//    std::vector<double> k_1_vec = {k_1, 0, 0, k_1};
//    std::vector<double> k_2_vec = {k_2*cos_phi_12, k_2*sin_phi_12, 0, k_2};

    std::vector<double> k_1_vec = {0,0,k_1,k_1};
    std::vector<double> k_2_vec = {k_2*sin_phi_12,0,k_2*cos_phi_12,k_2};

    std::vector<double> k_3_vec = q_ApB_vec(q_mA_vec(k_1_vec),q_mA_vec(k_2_vec));

    if (!is_vec_alright(k_1_vec))
    {
        std::cout << "NaN encountered in k_1_vec --> " << k_1_vec.at(0) << " " << k_1_vec.at(1) << " " << k_1_vec.at(2) << " " << k_1_vec.at(3) << std::endl;
    }

    if (!is_vec_alright(k_2_vec))
    {
        std::cout << "NaN encountered in k_2_vec --> " << k_2_vec.at(0) << " " << k_2_vec.at(1) << " " << k_2_vec.at(2) << " " << k_2_vec.at(3) << std::endl;
    }

    if (!is_vec_alright(k_3_vec))
    {
        std::cout << "NaN encountered in k_3_vec --> " << k_3_vec.at(0) << " " << k_3_vec.at(1) << " " << k_3_vec.at(2) << " " << k_3_vec.at(3) << std::endl;
    }

    //parameters in integrand
    params_B_1_loop_integrand args = {k_1_vec, k_2_vec, k_3_vec, z, class_obj, use_pk_nl};

    // in spherical coordinates: radius, mu (i.e. cos(theta) ), phi
    std::vector<double> lower_limits = { k_min_integral, -1, 0};
    std::vector<double> upper_limits = { k_max_integral,  1, 2*M_PI};

    hcubature_integration(B_1_loop_hcubature_integrand, static_cast<void *>(&args), lower_limits, upper_limits, 3, 4*1e4, result, error);

    return B_tree(k_1, k_2, k_3, z, class_obj, use_pk_nl) + result/(pow(2*M_PI,3));
}

// ######################################################################################

// matter bispectrum fitting formulae

double n(const double &k, const double &z, ClassEngine *class_obj)
{
    // is this precise enough?
    // symmetric difference quotient with h = 0.001 % --> the k term in both numerator and denominator are thus cancelled
    return 1 / class_obj->pk_lin(k, z) * (class_obj->pk_lin(k*1.001, z) - class_obj->pk_lin(k*0.999, z)) / (0.002);
}

double Q3(const double &n)
{
    double val = (4-pow(2,n)) / (1+pow(2,n+1));

    if (val < 0)
    {
        std::cout << "Q3 is negative with value = " << val << " for n = " << n << std::endl;
        return 0.0;
    }
    else
        return val;
}

double F2_eff(const double &k_1, const double &k_2, const double &k_3,
              const double &a_1, const double &a_2, const double &b_1, const double &b_2, const double &c_1, const double &c_2)
{
    if(k_1 == 0 || k_2 == 0)
        return 0.0;

    double cos_phi_12 = 0.5*(k_3*k_3-k_1*k_1-k_2*k_2)/(k_1*k_2); // law of cosines
    return 5.0/7.0*a_1*a_2 + 1/2.0*cos_phi_12*(k_1/k_2 + k_2/k_1)*b_1*b_2 + 2.0/7.0*pow(cos_phi_12,2)*c_1*c_2;
}

// ######################################################################################

// Scoccimarro and Couchman matter bispectrum fitting formula

double a_SC(const double &n, const double &q, const double &sigma8_z)
{
    return (1 + pow(sigma8_z,a6_SC)*sqrt(0.7*Q3(n))*pow(q*a1_SC,n+a2_SC)) / (1 + pow(q*a1_SC,n+a2_SC));
}

double b_SC(const double &n, const double &q)
{
    return (1 + 0.2*a3_SC*(n+3)*pow(q,n+3)) / (1 + pow(q,n+3.5));
}

double c_SC(const double &n, const double &q)
{
    return (1 + 4.5*a4_SC/(1.5 + pow(n+3,4))*pow(q*a5_SC,n+3)) / (1 + pow(q*a5_SC,n+3.5));
}

double B_SC(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    // Already ensure beforehand that the input k_i modes form a closed triangle

    double n_1 = class_obj->get_n_eff_from_lin_Pk(k_1, 0);
    double n_2 = class_obj->get_n_eff_from_lin_Pk(k_2, 0);
    double n_3 = class_obj->get_n_eff_from_lin_Pk(k_3, 0);

    double k_NL_z = class_obj->get_k_NL_from_lin_Pk_interp(z);

    double q_1 = k_1/k_NL_z;
    double q_2 = k_2/k_NL_z;
    double q_3 = k_3/k_NL_z;

    double sigma8_z = class_obj->get_sigma8_z(z);

    double Pk_1 = class_obj->pk(k_1,z,use_pk_nl);
    double Pk_2 = class_obj->pk(k_2,z,use_pk_nl);
    double Pk_3 = class_obj->pk(k_3,z,use_pk_nl);

    double a_1 = a_SC(n_1,q_1,sigma8_z);
    double a_2 = a_SC(n_2,q_2,sigma8_z);
    double a_3 = a_SC(n_3,q_3,sigma8_z);

    double b_1 = b_SC(n_1,q_1);
    double b_2 = b_SC(n_2,q_2);
    double b_3 = b_SC(n_3,q_3);

    double c_1 = c_SC(n_1,q_1);
    double c_2 = c_SC(n_2,q_2);
    double c_3 = c_SC(n_3,q_3);

    double B = 2.0*(F2_eff(k_1,k_2,k_3,a_1,a_2,b_1,b_2,c_1,c_2)*Pk_1*Pk_2
                   +F2_eff(k_2,k_3,k_1,a_2,a_3,b_2,b_3,c_2,c_3)*Pk_2*Pk_3
                   +F2_eff(k_3,k_1,k_2,a_3,a_1,b_3,b_1,c_3,c_1)*Pk_3*Pk_1);

    return B;
}

// ######################################################################################

// Gil-Marin matter bispectrum fitting formula

double a_GM(const double &n, const double &q, const double &sigma8_z)
{
    return (1 + pow(sigma8_z,a6_GM)*sqrt(0.7*Q3(n))*pow(q*a1_GM,n+a2_GM)) / (1 + pow(q*a1_GM,n+a2_GM));
}

double b_GM(const double &n, const double &q)
{
    return (1 + 0.2*a3_GM*(n+3)*pow(q*a7_GM,n+3+a8_GM)) / (1 + pow(q*a7_GM,n+3.5+a8_GM));
}

double c_GM(const double &n, const double &q)
{
    return (1 + 4.5*a4_GM/(1.5 + pow(n+3,4))*pow(q*a5_GM,n+3+a9_GM)) / (1 + pow(q*a5_GM,n+3.5+a9_GM));
}

double B_GM(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    // Already ensure beforehand that the input k_i modes form a closed triangle

//    double n_1 = n(k_1, z, class_obj);
//    double n_2 = n(k_2, z, class_obj);
//    double n_3 = n(k_3, z, class_obj);

    double n_1 = class_obj->get_n_eff_from_lin_Pk(k_1, 0);
    double n_2 = class_obj->get_n_eff_from_lin_Pk(k_2, 0);
    double n_3 = class_obj->get_n_eff_from_lin_Pk(k_3, 0);

    double k_NL_z = class_obj->get_k_NL_from_lin_Pk_interp(z);

    double q_1 = k_1/k_NL_z;
    double q_2 = k_2/k_NL_z;
    double q_3 = k_3/k_NL_z;

    //double sigma8_z = class_obj->get_sigma8_z(0); // potentially, this was used previously in the i3pt 2021 papers (both v1 and v2)
    double sigma8_z = class_obj->get_sigma8_z(z); // THIS IS CORRECT ---> CHANGE TO THIS! This doesn't make a huge difference compared to class_obj->get_sigma8_z(0)

    double Pk_1 = class_obj->pk(k_1,z,use_pk_nl);
    double Pk_2 = class_obj->pk(k_2,z,use_pk_nl);
    double Pk_3 = class_obj->pk(k_3,z,use_pk_nl);

    double a_1 = a_GM(n_1,q_1,sigma8_z);
    double a_2 = a_GM(n_2,q_2,sigma8_z);
    double a_3 = a_GM(n_3,q_3,sigma8_z);

    double b_1 = b_GM(n_1,q_1);
    double b_2 = b_GM(n_2,q_2);
    double b_3 = b_GM(n_3,q_3);

    double c_1 = c_GM(n_1,q_1);
    double c_2 = c_GM(n_2,q_2);
    double c_3 = c_GM(n_3,q_3);

    double B = 2.0*(F2_eff(k_1,k_2,k_3,a_1,a_2,b_1,b_2,c_1,c_2)*Pk_1*Pk_2
                   +F2_eff(k_2,k_3,k_1,a_2,a_3,b_2,b_3,c_2,c_3)*Pk_2*Pk_3
                   +F2_eff(k_3,k_1,k_2,a_3,a_1,b_3,b_1,c_3,c_1)*Pk_3*Pk_1);

    return B;
}

// ######################################################################################

// bihalofit fitting formula

double B_bihalofit(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    assert(compute_bihalofit);

    assert(class_obj->get_wa_fld() == 0);

    double D1 = class_obj->get_D_plus_z(z);
    double log10sigma8z = log10(D1*class_obj->get_sigma8_z(0));
    double ns = class_obj->get_n_s();
    double h = class_obj->get_h();

    double k1_h = k_1 / h; // [h/Mpc]
    double k2_h = k_2 / h;
    double k3_h = k_3 / h;

    double k[4], kt;
    k[1]=k1_h, k[2]=k2_h, k[3]=k3_h;
    for(int i=1;i<=3;i++) // sorting k[i] such that k[1]>=k[2]>=k[3]
    {
        for(int j=i+1;j<=3;j++)
        {
            if(k[i]<k[j])
            {
                kt=k[j];
                k[j]=k[i];
                k[i]=kt;
            }
        }
    }

    if(k[1]>k[2]+k[3])
    {
        //printf("Error: triangle is not closed! \n");
        return 0.0;
    }

    if( z > 10.)
        return B_tree(k_1,k_2,k_3,z,class_obj,false);

    double R_NL=class_obj->bihalofit_R_NL_interp(z); // Eq.(B1) 1st part R_NL
    double n_eff=class_obj->bihalofit_n_eff_interp(z); // Eq.(B2)

    double r1=k[3]/k[1];
    double r2=(k[2]+k[3]-k[1])/k[1];   // Eq.(B8)

    double q[4];
    q[1]=k1_h*R_NL, q[2]=k2_h*R_NL, q[3]=k3_h*R_NL;  // dimensionless wavenumbers

    // 1-halo term parameters in Eq.(B7)
    double an=pow(10.,-2.167-2.944*log10sigma8z-1.106*pow(log10sigma8z,2)-2.865*pow(log10sigma8z,3)-0.310*pow(r1,pow(10.,0.182+0.57*n_eff)));
    double bn=pow(10.,-3.428-2.681*log10sigma8z+1.624*pow(log10sigma8z,2)-0.095*pow(log10sigma8z,3));
    double cn=pow(10.,0.159-1.107*n_eff);
    double alphan=pow(10.,-4.348-3.006*n_eff-0.5745*pow(n_eff,2)+pow(10.,-0.9+0.2*n_eff)*pow(r2,2));
    if(alphan>1.-(2./3.)*ns)
        alphan=1.-(2./3.)*ns;
    double betan=pow(10.,-1.731-2.845*n_eff-1.4995*pow(n_eff,2)-0.2811*pow(n_eff,3)+0.007*r2);

    // 1-halo term bispectrum in Eq.(B4)
    double BS1h=1.;
    for(int i=1;i<=3;i++)
    {
        //if (q[i] == 0 )
        //    continue;
        BS1h*=1./(an*pow(q[i],alphan)+bn*pow(q[i],betan))/(1.+1./(cn*q[i]));
    }

    // 3-halo term parameters in Eq.(B9)
    double fn=pow(10.,-10.533-16.838*n_eff-9.3048*pow(n_eff,2)-1.8263*pow(n_eff,3));
    double gn=pow(10.,2.787+2.405*n_eff+0.4577*pow(n_eff,2));
    double hn=pow(10.,-1.118-0.394*n_eff);
    double mn=pow(10.,-2.605-2.434*log10sigma8z+5.71*pow(log10sigma8z,2));
    double nn=pow(10.,-4.468-3.08*log10sigma8z+1.035*pow(log10sigma8z,2));
    double mun=pow(10.,15.312+22.977*n_eff+10.9579*pow(n_eff,2)+1.6586*pow(n_eff,3));
    double nun=pow(10.,1.347+1.246*n_eff+0.4525*pow(n_eff,2));
    double pn=pow(10.,0.071-0.433*n_eff);
    double en=pow(10.,-0.632+0.646*n_eff);

    double PSE[4];
    for(int i=1;i<=3;i++)
    {
        //if (q[i] == 0 )
        //    continue;
        PSE[i]=(1.+fn*pow(q[i],2))/(1.+gn*q[i]+hn*pow(q[i],2))*pow(D1,2)*class_obj->bihalofit_pk_lin(q[i]/R_NL,0) +
                1./(mn*pow(q[i],mun)+nn*pow(q[i],nun))/(1.+pow(pn*q[i],-3));  // enhanced P(k,z) in Eq.(B6)
    }

    double dn=pow(10.,-0.483+0.892*log10sigma8z-0.086*class_obj->get_Omega_m_z(z));

    double F2_bihalofit_123 = F2_EdS(k1_h,k2_h,k3_h)+dn*q[3];
    double F2_bihalofit_231 = F2_EdS(k2_h,k3_h,k1_h)+dn*q[1];
    double F2_bihalofit_312 = F2_EdS(k3_h,k1_h,k2_h)+dn*q[2];

    // 3-halo term bispectrum in Eq.(B5)
    double BS3h = 2.0*( F2_bihalofit_123*PSE[1]*PSE[2] + F2_bihalofit_231*PSE[2]*PSE[3] + F2_bihalofit_312*PSE[3]*PSE[1] );
    for(int i=1;i<=3;i++)
    {
        //if (q[i] == 0 )
        //    continue;
        BS3h*=1./(1.+en*q[i]);
    }

    return (BS1h+BS3h)/pow(h,6);
}

double B_bihalofit_baryon_ratio(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj) // bihalofit bispectrum ratio with to without baryons
{
    assert(class_obj->get_wa_fld() == 0);

    double h = class_obj->get_h();

    double k1_h = k_1 / h; // [h/Mpc]
    double k2_h = k_2 / h;
    double k3_h = k_3 / h;

    if(z > 5.0)
        return 1.0;  // baryon_ratio is calibrated at z=0-5

    double a=1./(1.+z);

    double k[4],x[4];
    k[1]=k1_h, k[2]=k2_h, k[3]=k3_h;
    for(int i=1;i<=3;i++)
        x[i]=log10(k[i]);

    double A0;
    if(a > 0.5)
        A0=0.068*pow(a-0.5,0.47);
    else
        A0=0.;

    double mu0=0.018*a+0.837*a*a;
    double sigma0=0.881*mu0;
    double alpha0=2.346;

    double A1;
    if(a > 0.2)
        A1=1.052*pow(a-0.2,1.41);
    else
        A1=0.;

    double mu1=fabs(0.172+3.048*a-0.675*a*a);
    double sigma1=(0.494-0.039*a)*mu1;

    double ks=29.90-38.73*a+24.30*a*a;
    double alpha2=2.25;
    double beta2=0.563/(pow(a/0.060,0.02)+1.)/alpha2;

    double Rb=1.0;
    for(int i=1;i<=3;i++)
    {
        Rb*=A0*exp(-pow(fabs(x[i]-mu0)/sigma0,alpha0))-A1*exp(-pow(fabs(x[i]-mu1)/sigma1,2))+pow(1.+pow(k[i]/ks,alpha2),beta2);   // Eq.(C1)
    }

    return Rb;
}

// ######################################################################################

// response function squeezed matter bispectrum

double B_squeezed_RF(const double &k_h, const double &k_m, const double &k_s, const double &z, ClassEngine *class_obj, bool use_pk_nl)
{
    // Already ensure beforehand that the input k_i modes form a closed triangle AND are already sorted in descending order i.e. k_h > k_m > k_s

    // - negative sign in front to get cosine of the angle between the hard and soft vectors and not the triangular sides
    //double mu_h_s = - 0.5*(k_m*k_m + k_s*k_s - k_h*k_h) / (k_m*k_s); // cosine of the angle between k_m and k_s
    double mu_h_s = - 0.5*(k_h*k_h + k_s*k_s - k_m*k_m) / (k_h*k_s); // cosine of the angle between k_h and k_s - is this more correct?

    double n_h = class_obj->get_n_eff_from_nl_Pk(k_h,z);
    double R_1 = 1.0 - n_h/3.0 + class_obj->get_G_1_k_z_interp(k_h,z);
    double R_K = class_obj->get_G_K_k_z_interp(k_h,z) - n_h;

    double Pk_h = class_obj->pk(k_h,z,true);
//    double Pk_m = class_obj->pk(k_m,z,true);
    double Pk_s = class_obj->pk(k_s,z,true);

    double B = (R_1 + ( pow(mu_h_s,2) - 1.0/3.0 )*R_K )*Pk_h*Pk_s;

    return B;
}

// ######################################################################################
// ######################################################################################
// ######################################################################################

// primordial non-Gaussianity

// Local type

double B_primordial_local(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj)
{
    if (k_1 == 0 || k_2 == 0 || k_3 == 0)
        return 0.0;

    double M_k_1_z = class_obj->M_poisson_factor(k_1, z);
    double M_k_2_z = class_obj->M_poisson_factor(k_2, z);
    double M_k_3_z = class_obj->M_poisson_factor(k_3, z);

//    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
//    double Pk_1 = class_obj->pk_lin(k_1,0);
//    double Pk_2 = class_obj->pk_lin(k_2,0);
//    double Pk_3 = class_obj->pk_lin(k_3,0);

//    double B_primordial_local = 2*class_obj->get_f_NL_local()*(Pk_1*Pk_2*Mk_3_z/(Mk_1_z*Mk_2_z)
//                                                              +Pk_2*Pk_3*Mk_1_z/(Mk_2_z*Mk_3_z)
//                                                              +Pk_3*Pk_1*Mk_2_z/(Mk_3_z*Mk_1_z));

    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
    double P_phi_k_1 = class_obj->pk_gravitational_potential(k_1);
    double P_phi_k_2 = class_obj->pk_gravitational_potential(k_2);
    double P_phi_k_3 = class_obj->pk_gravitational_potential(k_3);

    double B_phi_primordial_local = 2*class_obj->get_f_NL_local()*(P_phi_k_1*P_phi_k_2 + P_phi_k_2*P_phi_k_3 + P_phi_k_3*P_phi_k_1);

    return B_phi_primordial_local*M_k_1_z*M_k_2_z*M_k_3_z;
}

double B_primordial_equilateral(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj)
{
    if (k_1 == 0 || k_2 == 0 || k_3 == 0)
        return 0.0;

    double M_k_1_z = class_obj->M_poisson_factor(k_1, z);
    double M_k_2_z = class_obj->M_poisson_factor(k_2, z);
    double M_k_3_z = class_obj->M_poisson_factor(k_3, z);

    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
    double P_phi_k_1 = class_obj->pk_gravitational_potential(k_1);
    double P_phi_k_2 = class_obj->pk_gravitational_potential(k_2);
    double P_phi_k_3 = class_obj->pk_gravitational_potential(k_3);

    double B_phi_primordial_equilateral = 6*class_obj->get_f_NL_equilateral()*(- (P_phi_k_1*P_phi_k_2 + P_phi_k_2*P_phi_k_3 + P_phi_k_3*P_phi_k_1)
                                                                               - 2*pow(P_phi_k_1*P_phi_k_2*P_phi_k_3, 2/3.0)
                                                                               + (pow(P_phi_k_1, 1/3.0)*pow(P_phi_k_2, 2/3.0)*P_phi_k_3 +
                                                                                  pow(P_phi_k_1, 1/3.0)*pow(P_phi_k_3, 2/3.0)*P_phi_k_2 +
                                                                                  pow(P_phi_k_2, 1/3.0)*pow(P_phi_k_1, 2/3.0)*P_phi_k_3 +
                                                                                  pow(P_phi_k_2, 1/3.0)*pow(P_phi_k_3, 2/3.0)*P_phi_k_1 +
                                                                                  pow(P_phi_k_3, 1/3.0)*pow(P_phi_k_1, 2/3.0)*P_phi_k_2 +
                                                                                  pow(P_phi_k_3, 1/3.0)*pow(P_phi_k_2, 2/3.0)*P_phi_k_1) );


    return B_phi_primordial_equilateral*M_k_1_z*M_k_2_z*M_k_3_z;
}

double B_primordial_orthogonal(const double &k_1, const double &k_2, const double &k_3, const double &z, ClassEngine *class_obj)
{
    if (k_1 == 0 || k_2 == 0 || k_3 == 0)
        return 0.0;

    double M_k_1_z = class_obj->M_poisson_factor(k_1, z);
    double M_k_2_z = class_obj->M_poisson_factor(k_2, z);
    double M_k_3_z = class_obj->M_poisson_factor(k_3, z);

    // Class object being used; k [1/Mpc] and P(k) [Mpc^3]
    double P_phi_k_1 = class_obj->pk_gravitational_potential(k_1);
    double P_phi_k_2 = class_obj->pk_gravitational_potential(k_2);
    double P_phi_k_3 = class_obj->pk_gravitational_potential(k_3);

    double B_phi_primordial_orthogonal = 6*class_obj->get_f_NL_orthogonal()*(- 3*(P_phi_k_1*P_phi_k_2 + P_phi_k_2*P_phi_k_3 + P_phi_k_3*P_phi_k_1)
                                                                             - 8*pow(P_phi_k_1*P_phi_k_2*P_phi_k_3, 2/3.0)
                                                                             + 3*(pow(P_phi_k_1, 1/3.0)*pow(P_phi_k_2, 2/3.0)*P_phi_k_3 +
                                                                                  pow(P_phi_k_1, 1/3.0)*pow(P_phi_k_3, 2/3.0)*P_phi_k_2 +
                                                                                  pow(P_phi_k_2, 1/3.0)*pow(P_phi_k_1, 2/3.0)*P_phi_k_3 +
                                                                                  pow(P_phi_k_2, 1/3.0)*pow(P_phi_k_3, 2/3.0)*P_phi_k_1 +
                                                                                  pow(P_phi_k_3, 1/3.0)*pow(P_phi_k_1, 2/3.0)*P_phi_k_2 +
                                                                                  pow(P_phi_k_3, 1/3.0)*pow(P_phi_k_2, 2/3.0)*P_phi_k_1) );


    return B_phi_primordial_orthogonal*M_k_1_z*M_k_2_z*M_k_3_z;
}
