#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <float.h>
#include <math.h>

#include "event_detection.h"

/**
 *   Compute cumulative sum and sum of squares for a vector of data
 *   data      double[d_length]   Data to be summed over (in)
 *   sum       double[d_length]   Vector to store sum (out)
 *   sumsq     double[d_length]   Vector to store sum of squares (out)
 *   d_length                     Length of data vector
 **/
void compute_sum_sumsq(const double *restrict data, double *restrict sum,
                       double *restrict sumsq, size_t d_length) {
    size_t i;

    // Basic contracts
    assert(NULL != data);
    assert(NULL != sum);
    assert(NULL != sumsq);
    assert(d_length > 0);

    sum[0] = data[0];
    sumsq[0] = data[0] * data[0];
    for (i = 1; i < d_length; ++i) {
        sum[i] = sum[i - 1] + data[i];
        sumsq[i] = sumsq[i - 1] + data[i] * data[i];
    }
}

/**
 *   Compute windowed t-statistic from summary information
 *   sum       double[d_length]  Cumulative sums of data (in)
 *   sumsq     double[d_length]  Cumulative sum of squares of data (in)
 *   tstat     double[d_length]  T-statistic (out)
 *   d_length                    Length of data vector
 *   w_length                    Window length to calculate t-statistic over
 **/
void compute_tstat(const double *restrict sum, const double *restrict sumsq,
                   double *restrict tstat, size_t d_length, size_t w_length,
                   bool pooled) {
    size_t i;
    const double eta = 1e-100;

    // Simple contracts
    assert(NULL != sum);
    assert(NULL != sumsq);
    assert(NULL != tstat);

    // Quick return:
    //   t-test not defined for number of points less than 2
    //   need at least as many points as twice the window length
    if (d_length < 2 * w_length || w_length < 2) {
        for (i = 0; i < d_length; ++i) {
            tstat[i] = 0.0;
        }
        return;
    }
    // fudge boundaries
    for (i = 0; i < w_length; ++i) {
        tstat[i] = 0;
        tstat[d_length - i - 1] = 0;
    }

    // get to work on the rest
    {
        double sum1, sum2, sumsq1, sumsq2, mean1, mean2, var1, var2;

        for (i = w_length; i <= d_length - w_length; ++i) {
            sum1 = sum[i - 1];
            sumsq1 = sumsq[i - 1];
            if (i > w_length) {
                sum1 -= sum[i - w_length - 1];
                sumsq1 -= sumsq[i - w_length - 1];
            }
            sum2 = sum[i + w_length - 1] - sum[i - 1];
            sumsq2 = sumsq[i + w_length - 1] - sumsq[i - 1];
            mean1 = sum1 / w_length;
            mean2 = sum2 / w_length;
            var1 = sumsq1 / w_length - mean1 * mean1;
            var2 = sumsq2 / w_length - mean2 * mean2;
            if (pooled) {
                var1 = (var1 + var2) / 2.0;
                var2 = var1;
            }
            // Prevent problem due to very small variances
            var1 = fmax(var1, eta);
            var2 = fmax(var2, eta);

            //t-stat
            //  Formula is a simplified version of Student's t-statistic for the
            //  special case where there are two samples of equal size with
            //  differing variance
            {
                const double delta = mean2 - mean1;
                const double totvar = var1 / w_length + var2 / w_length;
                tstat[i] = fabs(delta / sqrt(totvar));
            }
        }
    }
}

void short_long_peak_detector(DetectorPtr short_detector,
                              DetectorPtr long_detector,
                              const double peak_height, size_t * peaks) {
    size_t i, k;
    size_t peak_count = 0;
    DetectorPtr detector;
    DetectorPtr detectors[2] = { short_detector, long_detector };
    double current_value;

    assert(short_detector->signal_length == long_detector->signal_length);
    assert(NULL != peaks);

    for (i = 0; i < short_detector->signal_length; i++) {
        for (k = 0; k < 2; k++) {
            detector = detectors[k];
            //Carry on if we've been masked out
            if (detector->masked_to >= i) {
                continue;
            }

            current_value = detector->signal[i];

            if (detector->peak_pos == detector->DEF_PEAK_POS) {
                //CASE 1: We've not yet recorded a maximum
                if (current_value < detector->peak_value) {
                    //Either record a deeper minimum...
                    detector->peak_value = current_value;
                } else if (current_value - detector->peak_value > peak_height) {
                    // ...or we've seen a qualifying maximum
                    detector->peak_value = current_value;
                    detector->peak_pos = i;
                    //otherwise, wait to rise high enough to be considered a peak
                }
            } else {
                //CASE 2: In an existing peak, waiting to see if it is good
                if (current_value > detector->peak_value) {
                    //Update the peak
                    detector->peak_value = current_value;
                    detector->peak_pos = i;
                }
                //Dominate other tstat signals if we're going to fire at some point
                if (detector == short_detector) {
                    if (detector->peak_value > detector->threshold) {
                        long_detector->masked_to =
                            detector->peak_pos + detector->window_length;
                        long_detector->peak_pos = long_detector->DEF_PEAK_POS;
                        long_detector->peak_value = long_detector->DEF_PEAK_VAL;
                        long_detector->valid_peak = false;
                    }
                }
                //Have we convinced ourselves we've seen a peak
                if (detector->peak_value - current_value > peak_height
                    && detector->peak_value > detector->threshold) {
                    detector->valid_peak = true;
                }
                //Finally, check the distance if this is a good peak
                if (detector->valid_peak
                    && (i - detector->peak_pos) > detector->window_length / 2) {
                    //Emit the boundary and reset
                    peaks[peak_count] = detector->peak_pos;
                    peak_count++;
                    detector->peak_pos = detector->DEF_PEAK_POS;
                    detector->peak_value = current_value;
                    detector->valid_peak = false;
                }
            }
        }
    }
}

event_table detect_events(raw_table const rt)
{
    size_t const window1_length = 3;
    size_t const window2_length = 6;
    double const threshold1 = 1.4;
    double const threshold2 = 1.1;
    double const peak_height = 0.2;  // TODO(semen): pass on the cmd line

    double * data = calloc(rt.n, sizeof(double));
    double * sums = calloc(rt.n, sizeof(double));
    double * sumsqs = calloc(rt.n, sizeof(double));
    double * tstat1 = calloc(rt.n, sizeof(double));
    double * tstat2 = calloc(rt.n, sizeof(double));
    size_t * peaks = calloc(rt.n, sizeof(size_t));

    for( size_t i = 0; i < rt.n; ++i )
    {
        data[i] = (double) rt.raw[i];
    }

    compute_sum_sumsq(data, sums, sumsqs, rt.n);
    compute_tstat(sums, sumsqs, tstat1, rt.n, window1_length, false);
    compute_tstat(sums, sumsqs, tstat2, rt.n, window2_length, false);

    Detector short_detector = {
        .DEF_PEAK_POS = -1,
        .DEF_PEAK_VAL = 1e100,
        .signal = tstat1,
        .signal_length = rt.n,
        .threshold = threshold1,
        .window_length = window1_length,
        .masked_to = 0,
        .peak_pos = -1,
        .peak_value = 1e100,
        .valid_peak = false
    };

    Detector long_detector = {
        .DEF_PEAK_POS = -1,
        .DEF_PEAK_VAL = 1e100,
        .signal = tstat2,
        .signal_length = rt.n,
        .threshold = threshold2,
        .window_length = window2_length,
        .masked_to = 0,
        .peak_pos = -1,
        .peak_value = 1e100,
        .valid_peak = false
    };

    short_long_peak_detector(&short_detector, &long_detector, peak_height, peaks);

    size_t nevent = 0;
    for( size_t i = 0; i < rt.n; ++i )
    {
        if ( peaks[i] > 0 && peaks[i] < rt.n )
        {
            nevent++;
        }
    }
    printf("number of events = %ld\n", nevent);

    event_t * events = calloc(nevent, sizeof(event_t));

    double s = 0;
    double sm = 0;
    double smsq = 0;
    size_t j = 0;
    for( size_t i = 0; i < rt.n; ++i )
    {
        if ( peaks[i] > 0 && peaks[i] < rt.n )
        {
            size_t e = peaks[i];

            events[j].start = s;
            {
                double const ev_sample_len = (double) e - s;
                double const ev_mean = (sums[e-1] - sm) / ev_sample_len;
                double const variance = fmax( 0.0, (sumsqs[e-1] - smsq) / ev_sample_len - pow(ev_mean, 2));
                events[j].length = (float) ev_sample_len;
                events[j].mean = (float) ev_mean;
                events[j].stdv = (float) sqrt(variance);
            }
            events[j].pos = -1;
            events[j].state = -1;

            s = (double) e;
            sm = sums[e-1];
            smsq = sumsqs[e-1];

            j++;
        }
    }

    free(peaks);
    free(tstat2);
    free(tstat1);
    free(sumsqs);
    free(sums);
    free(data);

    return (event_table) {nevent, 0, nevent, events};
}
