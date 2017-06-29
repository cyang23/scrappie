#ifndef EVENT_DETECTION_H
#    define EVENT_DETECTION_H

#    include <stdlib.h>
#    include <stdint.h>

#    if defined(_MSC_VER)
#        define false   0
#        define true    1
#        define bool int
#        define _Bool int
#        define fmax max
#        define fmin min
#    else
#        include <stdbool.h>
#    endif

#    include "scrappie_structures.h"

typedef struct {
    int DEF_PEAK_POS;
    double DEF_PEAK_VAL;
    double *signal;
    size_t signal_length;
    double threshold;
    size_t window_length;
    size_t masked_to;
    int peak_pos;
    double peak_value;
    _Bool valid_peak;
} Detector;
typedef Detector *DetectorPtr;

void short_long_peak_detector(DetectorPtr short_detector,
                              DetectorPtr long_detector,
                              const double peak_height, size_t * peaks);

event_table detect_events(raw_table const rt);

#endif                          /* EVENT_DETECTION_H */
