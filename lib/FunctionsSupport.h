
#ifndef incl_FunctionsSupport_h
#define incl_FunctionsSupport_h


// some functions in support

extern "C" 
{ 
  int  decomplr_matrix_(double*, int*, int*, int*);

  void inverse_matrix_(double*, double*, int*, int*);

  void solve_matrix_(double*, int*, double*, double*, int*);

  void mult_matrix_p_(double*, double*, double*, int*, int*, int*);

  void unit_s_p_(double*);

  void set_s_f_(double*, double*, double*);
  void set_s_m_(double*, double*);
  void set_s_am_(double*, double*);

  void mult_u_ss_p_(double*, double*, double*);
  void mult_u_su_p_(double*, double*, double*);
  void mult_u_us_p_(double*, double*, double*);
  void mult_s_us_ap_(double*, double*, double*);
  void mult_s_su_ap_(double*, double*, double*);

  void mult_4_4s_f_(double*, double*, double*, double*);
  void set_4_p_(double*, double*);

  double norm2_s_(double*);

  void diff_tat_p_(double*, double*);

  double det_u_(double*);

  void pzero_(double*, int*);
 
  double dot_(double*, double*, int*);
}


#endif

