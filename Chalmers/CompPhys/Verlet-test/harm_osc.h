#ifndef HARM_OSC_H_ //A check so to not rerun all the definitions below.
#define HARM_OSC_H_

void stepper_euler2F(double m, double XV0[2], double XV[][2], double dt, int N_steps);
void stepper_verlet(double XV0[2], double XV[][2], double dt, int N_steps);

void save_to_file(char *filename, double XV[][2], double dt, int N_steps);

#endif
