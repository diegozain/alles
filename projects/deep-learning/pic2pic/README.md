# pic2pic: smooth to rough squares
diego domenzain
December 2020 @ Colorado School of Mines

__This is an example of mapping a smooth matrices to their rough equivalent.__

## Purpose

The purpose is to build __gradient boosters__ for a physics based inversion.

This project in itself is a proof of concept. 

## Map smooth squares to rough squares

_Output_

* __For training__. 24x24 matrix sub-divided in 3x3 blocks, each of size 8x8.
   * Each block has as value a number between 1/9 and 1.
* __In practice__. As the input, but with added high spatial-frequency content.

_Input_

* __For training__. Smooth version of the _output_ matrix (for training).
* __In practice__. A smooth image of something that should have more high spatial-frequency content.

## Problem set-up

1. Build 9! = 362880 matrices, each of size 24x24 and sub-divided in 3x3 blocks, each block of size 8x8. Each block has as value a number between 1/9 and 1.

1. Smooth each of the 9! matrices of step __1__.

1. The __learning machine__ will map the matrices of step __2__ to those of step __1__.

## Code

* The script ```idea_ex.m``` just plots the idea of what should happen. No training or model generation here.

* ```pic_generator.m``` builds the training data: ```x, y = b_, b```.
---

[![](../pics/pic2pic_idea.png)](./)

