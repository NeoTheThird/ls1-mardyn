/*
 * config.h
 *
 *  Created on: Jun 14, 2011
 *      Author: andreas
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#define MAX_BLOCK_SIZE 512

//#define REFERENCE_IMPLEMENTATION

// TODO: either use warpSize or WARP_SIZE!
#ifndef REFERENCE_IMPLEMENTATION
#define WARP_SIZE 32
#define NUM_WARPS 4
#else
#define WARP_SIZE 1
#define NUM_WARPS 1
#endif

#define BLOCK_SIZE (WARP_SIZE*NUM_WARPS)

#endif /* CONFIG_H_ */
