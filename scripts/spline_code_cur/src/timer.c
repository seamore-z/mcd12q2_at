/***************************************************************************//**
 * @file
 * @brief Functions for timing blocks of code
 ******************************************************************************/

#include <assert.h>   // assert()
#include <sys/time.h> // timeval, gettimeofday()

// local prototypes
float timeval_subtract (struct timeval, struct timeval);

/***************************************************************************//**
 * @brief Wallclock timer with subsecond resolution, starts, stops, and prints                                                                               
 * @param [in] start_timer Logical flag, start timer (1) or stop timer (0)
 * @return Elapsed time if start==1, else 0
 *
 * This function is used to start and stop a wallclock timer. To start, call
 * wallclock_timer(1). To stop and get elapsed time, call t=wallclock_timer(0).
 * The program will fail with an assertion error if you attempt to start or stop
 * more than once in a row.
 ******************************************************************************/
float wallclock_timer(int start_timer) {
   
    static int on;
    static struct timeval t0;
    struct timeval t1;

    if (start_timer == 1) {
        // start timer and exit
        assert(on == 0);
        assert(gettimeofday(&t0, 0) != -1);
        on = 1;
        return 0.0;

    } else {
        // stop timer and return elapsed wallclock time
        assert(on == 1);
        assert(gettimeofday(&t1, 0) != -1);
        on = 0;
        return timeval_subtract(t1, t0);
    }
}

/***************************************************************************//**
 * @brief Subtract two timevals, x-y, returning elapsed time in seconds
 * @param x [in] 
 * @param y [in]
 * @return Elapsed time in seconds
 * 
 * Modified from http://www.gnu.org/software/libc/manual/html_node/Date-and-Time.html
 ******************************************************************************/
float timeval_subtract(struct timeval x, struct timeval y) {

    int nsec, usec, sec;

    // Perform the carry for the later subtraction by updating y 
    if (x.tv_usec < y.tv_usec) {
        nsec = (y.tv_usec - x.tv_usec) / 1000000 + 1;
        y.tv_usec -= 1000000 * nsec;
        y.tv_sec += nsec;
    }
    if (x.tv_usec - y.tv_usec > 1000000) {
        nsec = (y.tv_usec - x.tv_usec) / 1000000;
        y.tv_usec += 1000000 * nsec;
        y.tv_sec -= nsec;
    }

    // compute the time difference
    sec = x.tv_sec - y.tv_sec;
    usec = x.tv_usec - y.tv_usec;

    return (float)sec + (float)usec/1000000.0;
}
