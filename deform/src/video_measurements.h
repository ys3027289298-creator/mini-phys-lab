/*
 * video_measurements.h
 *
 * Max squash ratios measured frame-by-frame off the art reference clip.
 *
 * *** PLACEHOLDER VALUES — REPLACE FROM THE REAL REFERENCE VIDEO ***
 *
 * The reference clip described in the task is not present in the workspace
 * (the only video on disk, ~/dock_take.mp4, is an unrelated UI docking
 * screen recording).  The numbers below are engineering placeholders that
 * follow the sqrt(drop-height) scaling of a linear spring contact, so the
 * full calibration pipeline can be built and exercised end to end.
 *
 * How to replace them once the clip arrives (60 fps assumed):
 *   1. ffmpeg -i ref.mp4 -vsync 0 frames/f%04d.png
 *   2. On the resting frame measure the vertical diameter D0 (pixels).
 *   3. Step through each impact, find the frame of maximum flattening,
 *      measure the vertical extent Dmin of the ball on that frame.
 *   4. epsilon = (D0 - Dmin) / D0, one value per drop height.
 *   5. Run tools/measure_video.py for an assisted pixel measurement,
 *      then paste the three values into VideoMeasure[] below.
 *
 * Only this file needs editing; make test and the calibration tool pick
 * the new values up automatically.
 */
#ifndef VIDEO_MEASUREMENTS_H
#define VIDEO_MEASUREMENTS_H

#define VIDEO_SOURCE          "PLACEHOLDER (reference clip not in workspace)"
#define VIDEO_FRAME_RATE_HZ   60
#define VIDEO_TOLERANCE       0.10f   /* allowed |sim - video| / video */

typedef struct {
    float drop_height;    /* m, floor to ball centre at release */
    float squash_ratio;   /* measured max squash, dimensionless  */
    int   impact_frame;   /* frame of maximum flattening, 0 = n/a */
} VideoMeasure;

/* Three drop heights.  Placeholder ratios scale as sqrt(h). */
static const VideoMeasure VideoMeasurements[] = {
    /* h [m]  epsilon   frame */
    {   0.5f,   0.080f,   0   },
    {   1.0f,   0.113f,   0   },
    {   2.0f,   0.160f,   0   },
};

#define VIDEO_N_MEASURES (int)(sizeof(VideoMeasurements) / sizeof(VideoMeasurements[0]))

#endif /* VIDEO_MEASUREMENTS_H */
