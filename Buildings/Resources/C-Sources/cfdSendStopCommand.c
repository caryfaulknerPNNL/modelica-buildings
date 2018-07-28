/*
 *
 * \file   cfdSendStopCommand.c
 *
 * \brief  Function to send a stop command to terminate the CFD simulation
 *
 * \author Wangda Zuo
 *         University of Miami
 *         W.Zuo@miami.edu
 *
 * \date   8/3/2013
 *
 */
#include "cfdCosimulation.h"

/*
 * Send a stop command to terminate the CFD simulation
 *
 * @return No return needed
 */
void cfdSendStopCommand(void *thread) {

size_t i = 0, imax = 10000;
char msg[100];

  cosim->para->flag = 0;
  ModelicaMessage("send stop command to FFD from destructor\n");

  /* Wait for the feedback from FFD*/
  while(cosim->para->flag==0 && i<imax) {
    if(cosim->para->ffdError==1) {
      ModelicaError(cosim->ffd->msg);
    }
    else {
      sleep(1);
	  sprintf(msg,"cosim->para->flag=%d and wait for %d s\n",cosim->para->flag,i);
	  ModelicaMessage(msg);
      i++;
    }
  }

  if(i<imax) {
    if(cosim->para->ffdError==1) {
      ModelicaError(cosim->ffd->msg);
    }
    else {
      ModelicaMessage("Successfully stopped the FFD simulation.\n");
    }
  }
  else {
    ModelicaError("Error: Cannot stop the FFD simulation in required time.");
  }

  free(cosim->para);
  free(cosim->modelica);
  free(cosim->ffd);
  free(cosim);
  
/* Windows*/
#ifdef _MSC_VER
	CloseHandle(thread);
/*  Linux*/
#else
  pthread_cancel(thread);
#endif

} /* End of cfdSendStopCommand*/
