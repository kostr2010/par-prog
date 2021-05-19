#ifndef CALC_H_INCLUDED
#define CALC_H_INCLUDED

#include <math.h>
#include <pthread.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

const size_t INTERVALS_N = 1e1;
const double EPSILON = 1e-3;

struct pair_t
{
    double first;
    double second;
};
typedef struct pair_t pair_t;

struct calc_space_t
{
    pair_t* intervals;
    size_t intervals_cur;

    pthread_mutex_t mux;
};
typedef struct calc_space_t calc_space_t;

double func(const double x)
{
    return sin(1.0 / x);
}

double func_second_der(const double x)
{
    const double xx = x * x;

    return (2 * x * cos(1 / x) - sin(1 / x)) / (xx * xx);
}

// only suitable for sin(1/x) function
double FindMaximum(const pair_t interval, const double epsilon)
{
    double max = fabs(func_second_der(interval.first));
    for (double i = interval.first; i < interval.second; i += 1e-6)
    {
        double res = fabs(func_second_der(i));
        if (res >= max)
        {
            max = res;
        }
        else
        {
            return max;
        }
    }

    return max;
}

double GetStep(const pair_t interval, const double epsilon)
{
    const double length = interval.second - interval.first;

    /**
     * lazy adaptive step
     */
    static const double DIV = 1e5;
    const double interm = fabs(interval.second + interval.first) / 2.0;
    const double res = interm - 1.0 / (PI + 1.0 / interm);

    if (res <= 2 * DIV * __DBL_EPSILON__)
    {
        return 2 * __DBL_EPSILON__;
    }
    else
    {
        return (res > length) ? (length / DIV) : (res / DIV);
    }

    // double step_suggested =
    //     sqrt(epsilon * 12 / ((interval.second - interval.first) * FindMaximum(interval,
    //     epsilon)));

    // return (step_suggested > length) ? (length) : (step_suggested);
}

double CalculatePartialIntegral(const pair_t interval, const double epsilon)
{
    const double step = GetStep(interval, epsilon);
    double sum = 0.0;

    for (size_t cur = 1; interval.first + cur * step < interval.second - step; cur++)
    {
        sum += func(interval.first + cur * step);
    }

    sum += func(interval.first) / 2.0 + func(interval.second - step) / 2.0;
    sum *= step;

    // printf("%f %f %.16f, max: %f ,%f\n",
    //        interval.first,
    //        interval.second,
    //        step,
    //        FindMaximum(interval, epsilon),
    //        sum);

    return sum;
}

void* Routine(void* arg)
{
    double* res = (double*)calloc(sizeof(double), 1);

    pthread_mutex_lock(&(((calc_space_t*)arg)->mux));
    size_t intervals_cur = ((calc_space_t*)arg)->intervals_cur++;

    if (intervals_cur >= INTERVALS_N)
    {
        pthread_mutex_unlock(&(((calc_space_t*)arg)->mux));

        return res;
    }

    pair_t* intervals = ((calc_space_t*)arg)->intervals;
    pthread_mutex_unlock(&(((calc_space_t*)arg)->mux));

    while (intervals_cur < INTERVALS_N)
    {
        const pair_t interval = intervals[intervals_cur];

        *res += CalculatePartialIntegral(interval, EPSILON);

        pthread_mutex_lock(&(((calc_space_t*)arg)->mux));
        intervals_cur = ((calc_space_t*)arg)->intervals_cur++;
        pthread_mutex_unlock(&(((calc_space_t*)arg)->mux));
    }

    return res;
}

#endif