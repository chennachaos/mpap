
#ifndef incl_SolverTime_h
#define incl_SolverTime_h



class SolverTime
{
  public:

    SolverTime(void);

    ~SolverTime();

    double factorise, solve, factoriseAndSolve, total;

    void print(void);

    void reset(void);
};







#endif



