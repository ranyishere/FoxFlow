
//simulation time in seconds
settings.TOTAL_TIME = 20.0;//25.0;//20.0;
settings.DELTA = 0.5/8.0; //unit of seconds

//The internal step of the solver should be at least smaller than delta
settings.DELTA_DELTA_T = settings.DELTA / 20.0;
settings.DELTA_T_MIN = settings.DELTA_DELTA_T;
settings.NUM_STEPS = settings.TOTAL_TIME / settings.DELTA;
