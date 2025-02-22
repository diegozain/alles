#include <stdio.h>
#include <unistd.h>  // for sleep() function

// gcc -o animac animac.c
// https://stackoverflow.com/questions/4842424/list-of-ansi-color-escape-sequences
// printf("\033[H\033[J"); // clear the screen

int main() {
  while (1) {
    // frame 1
    printf("\033[31m"); // 🔴 font
    printf("\n          !\n");
    printf("       \\(o_o)\n");
    printf("       _/   \\_\n");
    printf("\033[0m"); // normal color font
    fflush(stdout);
    sleep(1); // seconds

    printf("\r");     // move to begining of line
    printf("\033[A"); // go up one line
    printf("\033[A"); // go up one line
    printf("\033[A"); // go up one line
    printf("\033[A"); // go up one line

    // frame 2
    printf("\033[33m"); // 🟡 font
    // printf("\033[38;5;206;48;5;57m"); // some weird shit tbh
    printf("\n          .\n");
    printf("       _(._.)\n");
    printf("       _/   \\_\n");
    fflush(stdout);
    sleep(1); // seconds

    printf("\r");     // move to begining of line
    printf("\033[A"); // go up one line
    printf("\033[A"); // go up one line
    printf("\033[A"); // go up one line
    printf("\033[A"); // go up one line
  }

  return 0;
}
