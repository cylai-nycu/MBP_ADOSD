
double exp_fast_lb(double a);
double exp_fast_ub(double a);
double exp_fast(double a);
double better_exp_fast(double a);
double exp_fast_schraudolph(double a);


float expf_fast_ub(float a);
float expf_fast(float a);
double better_expf_fast(float a);
float expf_fast_lb(float a);

double log_fast_ankerl(double a);
double log_fast_ub(double a);
double log_fast(double a);
double log_fast_lb(double a);

float logf_fast_ub(float a);
float logf_fast(float a);
float logf_fast_lb(float a);

double pow_fast_ankerl(double a, double b);

float powf_fast(float a, float b);
float powf_fast_lb(float a, float b);
float powf_fast_ub(float a, float b);

double pow_fast_ub(double a, double b);
double pow_fast(double a, double b);
double pow_fast_lb(double a, double b);
double pow_fast_precise_ankerl(double a, double b);
double pow_fast_precise(double a, double b);
double better_pow_fast_precise(double a, double b);

float powf_fast_precise(float a, float b);
float better_powf_fast_precise(float a, float b);

void fast_math_test(void);
