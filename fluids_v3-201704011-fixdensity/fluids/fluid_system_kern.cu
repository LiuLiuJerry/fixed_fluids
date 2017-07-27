/*
  FLUIDS v.3 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2012-2013. Rama Hoetzlein, http://fluids3.com

  Attribute-ZLib license (* See additional part 4)

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  4. Any published work based on this code must include public acknowledgement
     of the origin. This includes following when applicable:
	   - Journal/Paper publications. Credited by reference to work in text & citation.
	   - Public presentations. Credited in at least one slide.
	   - Distributed Games/Apps. Credited as single line in game or app credit page.	 
	 Retaining this additional license term is required in derivative works.
	 Acknowledgement may be provided as:
	   Publication version:  
	      2012-2013, Hoetzlein, Rama C. Fluids v.3 - A Large-Scale, Open Source
	 	  Fluid Simulator. Published online at: http://fluids3.com
	   Single line (slides or app credits):
	      GPU Fluids: Rama C. Hoetzlein (Fluids v3 2013)

 Notes on Clause 4:
  The intent of this clause is public attribution for this contribution, not code use restriction. 
  Both commerical and open source projects may redistribute and reuse without code release.
  However, clause #1 of ZLib indicates that "you must not claim that you wrote the original software". 
  Clause #4 makes this more specific by requiring public acknowledgement to be extended to 
  derivative licenses. 

*/

#define CUDA_KERNEL
#include "fluid_system_kern.cuh"

#include "cutil_math.h"

#include "radixsort.cu"						// Build in RadixSort

__constant__ FluidParams		simData;
__constant__ uint				gridActive;

__global__ void insertParticles ( bufList buf, int pnum )
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;

	register float3 gridMin = simData.gridMin; //方盒的最低点
	register float3 gridDelta = simData.gridDelta; //方盒的大小
	register int3 gridRes = simData.gridRes;
	register int3 gridScan = simData.gridScanMax; //扫描的最大值
	register float poff = simData.psmoothradius / simData.psimscale;//********???

	register int		gs;
	register float3		gcf;
	register int3		gc;

	gcf = (buf.mpos[i] - gridMin) * gridDelta; //方盒起点到目前粒子位置的向量
	gc = make_int3( int(gcf.x), int(gcf.y), int(gcf.z) ); //变成整数
	gs = (gc.y * gridRes.z + gc.z)*gridRes.x + gc.x; //在整个网格中的序号
	if ( gc.x >= 1 && gc.x <= gridScan.x && gc.y >= 1 && gc.y <= gridScan.y && gc.z >= 1 && gc.z <= gridScan.z ) {
		buf.mgcell[i] = gs;											// Grid cell insert.把粒子数据映射到二维blocks中
		buf.mgndx[i] = atomicAdd ( &buf.mgridcnt[ gs ], 1 );		// Grid counts.// mgridcnt[gs]+1，mgndx记录mgridcnt原来的值

		gcf = (-make_float3(poff,poff,poff) + buf.mpos[i] - gridMin) * gridDelta;
		gc = make_int3( int(gcf.x), int(gcf.y), int(gcf.z) );
		gs = ( gc.y * gridRes.z + gc.z)*gridRes.x + gc.x;		
	} else {
		buf.mgcell[i] = GRID_UNDEF;		
	}
}

// the mutex variable
__device__ int g_mutex = 0;

// GPU simple synchronization function
__device__ void __gpu_sync(int goalVal)
{

	__threadfence ();

	// only thread 0 is used for synchronization
	if (threadIdx.x == 0) 
		atomicAdd(&g_mutex, 1);
	
	// only when all blocks add 1 to g_mutex will
	// g_mutex equal to goalVal
	while(g_mutex < goalVal) {			// infinite loop until g_mutx = goalVal
	}

	if ( blockIdx.x == 0 && threadIdx.x == 0 ) g_mutex = 0;
	
	__syncthreads();
}

// countingSortInPlace -- GPU_SYNC DOES NOT WORK
/*uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;		// particle index				
	if ( i >= pnum ) { __gpu_sync ( 2 ); return; }

	register float3	ipos, ivel, iveleval, iforce;
	register float	ipress, idens;
	register int	icell, indx, iclr;

	icell = buf.mgcell [ i ];
	indx = buf.mgndx [ i ];
	int sort_ndx = buf.mgridoff[ icell ] + indx;				// global_ndx = grid_cell_offet + particle_offset
	if ( icell == GRID_UNDEF ) { __gpu_sync ( 2 ); return; }

	ipos = buf.mpos [ i ];
	ivel = buf.mvel [ i ];
	iveleval = buf.mveleval [ i ];
	iforce = buf.mforce [ i ];
	ipress = buf.mpress [ i ];
	idens = buf.mdensity [ i ];
	iclr = buf.mclr [ i ];

	__gpu_sync ( 2 ) ; //threadfence();			// make sure every thread in all blocks has their data

	
	buf.mpos [ sort_ndx ] = ipos;
	buf.mvel [ sort_ndx ] = ivel;
	buf.mveleval [ sort_ndx ] = iveleval;
	buf.mforce [ sort_ndx ] = iforce;
	buf.mpress [ sort_ndx ] = ipress;
	buf.mdensity [ sort_ndx ] = idens;
	buf.mclr [ sort_ndx ] = iclr;*/



// Counting Sort - Index
__global__ void countingSortIndex ( bufList buf, int pnum )
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;		// particle index				
	if ( i >= pnum ) return;

	uint icell = buf.mgcell[i];
	uint indx =  buf.mgndx[i];
	int sort_ndx = buf.mgridoff[ icell ] + indx;				// global_ndx = grid_cell_offet + particle_offset
	if ( icell != GRID_UNDEF ) {
		buf.mgrid[ sort_ndx ] = i;					// index sort, grid refers to original particle order
	}
}

// Counting Sort - Full (deep copy)
__global__ void countingSortFull ( bufList buf, int pnum )
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;		// particle index	mul24忽略整数的最高八位，因为这里的整数都不太大，应该是加快速度吧			
	if ( i >= pnum ) return; //防止越界

	// Copy particle from original, unsorted buffer (msortbuf),
	// into sorted memory location on device (mpos/mvel)
	uint icell = *(uint*) (buf.msortbuf + pnum*BUF_GCELL + i*sizeof(uint) );
	uint indx =  *(uint*) (buf.msortbuf + pnum*BUF_GNDX + i*sizeof(uint) );		

	if ( icell != GRID_UNDEF ) {	  
		// Determine the sort_ndx, location of the particle after sort
	    int sort_ndx = buf.mgridoff[ icell ] + indx;				// global_ndx = grid_cell_offet + particle_offset	
		
		// Find the original particle data, offset into unsorted buffer (msortbuf)
		char* bpos = buf.msortbuf + i*sizeof(float3);

		// Transfer data to sort location
		buf.mgrid[ sort_ndx ] = sort_ndx;			// full sort, grid indexing becomes identity		
		buf.mpos[ sort_ndx ] =		*(float3*) (bpos);
		buf.mvel[ sort_ndx ] =		*(float3*) (bpos + pnum*BUF_VEL );
		buf.mveleval[ sort_ndx ] =	*(float3*) (bpos + pnum*BUF_VELEVAL );
		buf.mforce[ sort_ndx ] =	*(float3*) (bpos + pnum*BUF_FORCE );
		buf.mpress[ sort_ndx ] =	*(float*) (buf.msortbuf + pnum*BUF_PRESS + i*sizeof(float) );
		buf.mdensity[ sort_ndx ] =	*(float*) (buf.msortbuf + pnum*BUF_DENS + i*sizeof(float) );
		buf.mclr[ sort_ndx ] =		*(uint*) (buf.msortbuf + pnum*BUF_CLR+ i*sizeof(uint) );		// ((uint) 255)<<24; -- dark matter
		buf.mgcell[ sort_ndx ] =	icell;
		buf.mgndx[ sort_ndx ] =		indx;		
	}
}

// ***** UNUSED CODE (not working) ******
__global__ void countActiveCells ( bufList buf, int pnum )
{	
	if ( threadIdx.x == 0 ) {		
		// use only one processor
		
		//gridActive = -1;

		int last_ndx = buf.mgridoff [ simData.gridTotal-1 ] + buf.mgridcnt[ simData.gridTotal-1 ] - 1;
		int last_p = buf.mgrid[ last_ndx ];
		int last_cell = buf.mgcell[ last_p ];
		int first_p = buf.mgrid[ 0 ];
		int first_cell = buf.mgcell[ first_p ] ;

		int cell, cnt = 0, curr = 0;
		cell = first_cell;
		while ( cell < last_cell ) {			
			buf.mgridactive[ cnt ] = cell;			// add cell to active list
			cnt++;
			curr += buf.mgridcnt[cell];				// advance to next active cell
			// id = buf.mgrid[curr];				// get particle id -- when unsorted only
			cell = buf.mgcell [ curr ];				// get cell we are in -- use id when unsorted
		}
		// gridActive = cnt;
	}
	__syncthreads();
}


__device__ float contributePressure ( int i, float3 p, int cell, bufList buf )
{			
	float3 dist;
	float dsq, c, sum;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2 / d2;
	
	sum = 0.0;

	if ( buf.mgridcnt[cell] == 0 ) return 0.0;
	
	int cfirst = buf.mgridoff[ cell ];
	int clast = cfirst + buf.mgridcnt[ cell ];
	
	for ( int cndx = cfirst; cndx < clast; cndx++ ) {
		dist = p - buf.mpos[ buf.mgrid[cndx] ];
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		if ( dsq < r2 && dsq > 0.0) {  //判断粒子的距离
			c = (r2 - dsq)*d2;
			sum += c * c * c;				
		} 
	}
	
	return sum;
}
			
__global__ void computePressure ( bufList buf, int pnum )
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;

	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[ i ];  //第i个粒子所在的网格
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

	// Sum Pressures
	float3 pos = buf.mpos[ i ];
	float sum = 0.0;
	for (int c=0; c < simData.gridAdjCnt; c++) {//周围的gridAdjCnt个网格
		sum += contributePressure ( i, pos, gc + simData.gridAdj[c], buf );//第c个网格的地址
	}
	__syncthreads();
		
	// Compute Density & Pressure
	sum = sum * simData.pmass * simData.poly6kern; //此处得到的是密度
	if ( sum == 0.0 ) sum = 1.0;
	buf.mpress[ i ] = ( sum - simData.prest_dens ) * simData.pintstiff;
	buf.mdensity[ i ] = 1.0f / sum;
}

		
__global__ void computeQuery ( bufList buf, int pnum )
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;

	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[ i ];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

	// Sum Pressures
	float3 pos = buf.mpos[ i ];
	float sum = 0.0;
	for (int c=0; c < simData.gridAdjCnt; c++) {
		sum += 1.0;
	}
	__syncthreads();
	
}

/*FindNeighbors
int cid = blockIdx.x * blockSize.x + blockIdx.y;   // cluster id	
int pid = threadIdx.x;		           // 0 to 85 (max particles per cell)	
__shared__ Particle  clist[ 85 ];	
__shared__ Particle  plist[ 85*8 ];
if ( pid < clusterCnt[cid] )  
	clist [ pid ] = particles [ clusterNdx[cid] + pid ];

for ( gid = 0;  gid < 8;  gid++ ) {
	if ( pid < gridCnt[  cid + group[gid] ] )  
		plist [ cid*CELL_CNT + pid ] = particles [ sortNdx[ cid + group[gid] ]  + pid ]; 	}

__syncthreads();	
	
for ( int j = 0; j < cellcnt;  j++ ) {
	dst = plist[ pid ] - plist[ j ];
	if ( dst < R2 ) {
     		  ...
	}
}*/

/*grid		    block
<gx, gy, gz>    <1, 32, 64>
256, 256, 256  
total:  */


#define LOCAL_PMAX		896
#define NUM_CELL		27
#define LAST_CELL		26
#define CENTER_CELL		13

__global__ void computePressureGroup ( bufList buf, int pnum )
{
	__shared__ float3	cpos[ LOCAL_PMAX ];

	__shared__ int		ncnt[ NUM_CELL ];
	__shared__ int		ngridoff[ NUM_CELL ];
	__shared__ int		noff[ NUM_CELL ];
	
	int bid = __mul24( blockIdx.y, gridDim.x ) + blockIdx.x;
	if ( bid > gridActive ) return;				// block must be in a valid grid
	uint cell = buf.mgridactive [ bid ];		// get grid cell (from blockID 1:1)
	register int i = -1;
	register float3 ipos;

	uint ndx = threadIdx.x;							
	if ( ndx < buf.mgridcnt[cell] ) {
		i = buf.mgridoff[cell] + ndx;		// particle id to process
		ipos = buf.mpos[ i ];
	}
	int gid = threadIdx.x;

	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2 / d2;
	register float3 dist;
	register float c, dsq, sum;
	int neighbor;

	// copy neighbor cell counts to shared mem
	if ( gid < NUM_CELL ) {
		int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
		neighbor = cell - nadj + simData.gridAdj[gid];					// neighbor cell id
		ncnt[gid] = buf.mgridcnt [ neighbor ];	
		ngridoff[gid] = buf.mgridoff [ neighbor ];
	}
	__syncthreads ();

	if ( gid == 0 ) {									// compute neighbor local ndx (as prefix sum)
		int nsum = 0;
		for (int z=0; z < NUM_CELL; z++) {				// 27-step prefix sum
			noff[z] = nsum;
			nsum += ncnt[z];
		}
	}
	__syncthreads ();

	// copy particles into shared memory
	if ( gid < NUM_CELL ) {
		for (int j=0; j < ncnt[gid]; j++ ) {
			neighbor = buf.mgrid [ ngridoff[gid] + j ];		// neighbor particle id
			ndx = noff[ gid ] + j;
			cpos[ ndx ] = buf.mpos [ neighbor ];
		}
	}
	__syncthreads ();

	
	// compute pressure for current particle
	if ( i == -1 ) return;
	
	int jnum = noff[LAST_CELL] + ncnt[LAST_CELL];
	sum = 0.0;
	for (int j = 0; j < jnum; j++) {
		dist = ipos - cpos[ j ];
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);			
		if ( dsq > 0.0 && dsq < r2 ) {
			c = (r2 - dsq)*d2;
			sum += c * c * c;
		}
	}	
	__syncthreads ();

	// put result into global mem
	sum = sum * simData.pmass * simData.poly6kern;
	if ( sum == 0.0 ) sum = 1.0;
	buf.mpress[ i ] = ( sum - simData.prest_dens ) * simData.pintstiff;
	buf.mdensity[ i ] = 1.0f / sum; 	
}


__device__ float3 contributeForce ( int i, float3 ipos, float3 iveleval, float ipress, float idens, int cell, bufList buf )
{			
	float dsq, c;	
	float pterm;
	float3 dist, force;	
	int j;					

	if ( buf.mgridcnt[cell] == 0 ) return make_float3(0,0,0);	

	force = make_float3(0,0,0);

	for ( int cndx = buf.mgridoff[ cell ]; cndx < buf.mgridoff[ cell ] + buf.mgridcnt[ cell ]; cndx++ ) {										
		j = buf.mgrid[ cndx ];				
		dist = ( ipos - buf.mpos[ j ] );		// dist in cm   两点间距离
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		if ( dsq < simData.rd2 && dsq > 0) {			
			dsq = sqrt(dsq * simData.d2);
			c = ( simData.psmoothradius - dsq ); 
			pterm = simData.psimscale * -0.5f * c * simData.spikykern * ( ipress + buf.mpress[ j ] ) / dsq;		
			//就不乘以m了，反正算加速度还要除回来...但是乘以idens是干嘛呀，既然是倒数的话
			force += ( pterm * dist + simData.vterm * ( buf.mveleval[ j ] - iveleval )) * c * idens * (buf.mdensity[ j ] );
		}	
	}
	return force;
}


__global__ void computeForce ( bufList buf, int pnum)
{			
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;

	// Get search cell	
	uint gc = buf.mgcell[ i ];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;

	// Sum Pressures	
	register float3 force;
	force = make_float3(0,0,0);		

	for (int c=0; c < simData.gridAdjCnt; c++) {
		force += contributeForce ( i, buf.mpos[ i ], buf.mveleval[ i ], buf.mpress[ i ], buf.mdensity[ i ], gc + simData.gridAdj[c], buf );
	}
	buf.mforce[ i ] = force;
}
	

/*__global__ void computeForceNbr ( char* bufPnts, int* bufGrid, int numPnt )
{		
	uint ndx = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index		
	if ( ndx >= numPnt ) return;
				
	char* ioffs = bufPnts + __mul24(ndx, simData.stride );
	float3 ipos = *(float3*)	(ioffs + OFFSET_POS);
	float3 ivelval = *(float3*)	(ioffs + OFFSET_VELEVAL);
	float press = *(float*)		(ioffs + OFFSET_PRESS);
	float dens =  *(float*)		(ioffs + OFFSET_DENS);
	int icnt =  *(int*)			(ioffs + OFFSET_NBRCNT);

	char* joffs;
	float3 jpos, jveleval;

	float3 dist, force;		
	float c, ndistj, pterm, dterm, vterm;
		
	vterm = simData.lapkern * simData.visc;
		
	force = make_float3(0,0,0);
	for (int nbr=0; nbr < icnt; nbr++) {		// base 1, n[0] = count
		ndistj = bufNdist[ndx][nbr];
		joffs = bufPnts + __mul24(bufNeighbor[ndx][nbr], simData.stride);
		jpos = *(float3*)		(joffs + OFFSET_POS);
		jveleval = *(float3*)	(joffs + OFFSET_VELEVAL);
		c = ( simData.smooth_rad - ndistj ); 
		dist.x = ( ipos.x - jpos.x );		// dist in cm
		dist.y = ( ipos.y - jpos.y );
		dist.z = ( ipos.z - jpos.z );			
		pterm = simData.sim_scale * -0.5f * c * simData.spikykern * ( press + *(float*)(joffs+OFFSET_PRESS) ) / ndistj;
		dterm = c * dens * *(float*)(joffs+OFFSET_DENS);	
		force.x += ( pterm * dist.x + vterm * ( jveleval.x - ivelval.x )) * dterm;
		force.y += ( pterm * dist.y + vterm * ( jveleval.y - ivelval.y )) * dterm;
		force.z += ( pterm * dist.z + vterm * ( jveleval.z - ivelval.z )) * dterm;			
	}
	*(float3*) ( ioffs + OFFSET_FORCE ) = force;		
}*/

		
__global__ void advanceParticles ( float time, float dt, float ss, bufList buf, int numPnts )
{		
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if ( i >= numPnts ) return;
	
	if ( buf.mgcell[i] == GRID_UNDEF ) {
		buf.mpos[i] = make_float3(-1000,-1000,-1000);
		buf.mvel[i] = make_float3(0,0,0);
		return;
	}
			
	// Get particle vars
	register float3 accel, norm;
	register float diff, adj, speed;
	register float3 pos = buf.mpos[i];
	register float3 veval = buf.mveleval[i];

	// Leapfrog integration						
	accel = buf.mforce[i]; //只有这个force？
	accel *= simData.pmass; //本来公式中应该乘以两个mass，之前没乘，这次乘了一个，就不用除回去了
		
	// Boundaries
	// Y-axis
	//由于有一定的倾斜，所以要拿最低边界加上倾斜的那一点
	diff = simData.pradius - (pos.y - (simData.pboundmin.y + (pos.x-simData.pboundmin.x)*simData.pground_slope )) * ss;
	if ( diff > EPSILON ) {
		norm = make_float3( -simData.pground_slope, 1.0 - simData.pground_slope, 0);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval );
		norm *= adj; accel += norm;
	}

	diff = simData.pradius - ( simData.pboundmax.y - pos.y )*ss;
	if ( diff > EPSILON ) {
		norm = make_float3(0, -1, 0); //受到一个反弹的力
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval ); //***************？？？
		norm *= adj; accel += norm;//调整加速度
	}

	// X-axis    力的大小成正弦函数
	diff = simData.pradius - (pos.x - (simData.pboundmin.x + (sin(time*simData.pforce_freq)+1)*0.5 * simData.pforce_min))*ss; //X轴以前受到一个力，产生波浪
	if ( diff > EPSILON ) {
		norm = make_float3( 1, 0, 0);
		adj = (simData.pforce_min+1) * simData.pextstiff * diff - simData.pdamp * dot(norm, veval );
		norm *= adj; accel += norm;
	}
	diff = simData.pradius - ( (simData.pboundmax.x - (sin(time*simData.pforce_freq)+1)*0.5*simData.pforce_max) - pos.x)*ss;
	if ( diff > EPSILON ) {
		norm = make_float3(-1, 0, 0);
		adj = (simData.pforce_max+1) * simData.pextstiff * diff - simData.pdamp * dot(norm, veval );
		norm *= adj; accel += norm;
	}

	// Z-axis
	diff = simData.pradius - (pos.z - simData.pboundmin.z ) * ss; //z轴的上下边界
	if ( diff > EPSILON ) {
		norm = make_float3( 0, 0, 1 );
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval );
		norm *= adj; accel += norm;
	}
	diff = simData.pradius - ( simData.pboundmax.z - pos.z )*ss;
	if ( diff > EPSILON ) {
		norm = make_float3( 0, 0, -1 );
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval );
		norm *= adj; accel += norm;
	}
		
	// Gravity******************************************
	accel += simData.pgravity;

	// Accel Limit
	speed = accel.x*accel.x + accel.y*accel.y + accel.z*accel.z;
	if ( speed > simData.AL2 ) {
		accel *= simData.AL / sqrt(speed);
	}

	// Velocity Limit
	float3 vel = buf.mvel[i];
	speed = vel.x*vel.x + vel.y*vel.y + vel.z*vel.z;
	if ( speed > simData.VL2 ) {
		speed = simData.VL2;
		vel *= simData.VL / sqrt(speed);
	}

	// Ocean colors
	if ( speed > simData.VL2*0.2) {
		adj = simData.VL2*0.2;
		buf.mclr[i] += ((  buf.mclr[i] & 0xFF) < 0xFD ) ? +0x00000002 : 0;		// decrement R by one
		buf.mclr[i] += (( (buf.mclr[i]>>8) & 0xFF) < 0xFD ) ? +0x00000200 : 0;	// decrement G by one
		buf.mclr[i] += (( (buf.mclr[i]>>16) & 0xFF) < 0xFD ) ? +0x00020000 : 0;	// decrement G by one
	}
	if ( speed < 0.03 ) {		
		int v = int(speed/.01)+1;
		buf.mclr[i] += ((  buf.mclr[i] & 0xFF) > 0x80 ) ? -0x00000001 * v : 0;		// decrement R by one
		buf.mclr[i] += (( (buf.mclr[i]>>8) & 0xFF) > 0x80 ) ? -0x00000100 * v : 0;	// decrement G by one
	}
	////My Add*************************************************************************************************
	register int phase = buf.mphase[i];
	if (phase < 1){
		buf.mclr[i] = 0xFBFBFBFB;
	}
	////*******************************************************************************************************
	//
	////-- surface particle density 
	//buf.mclr[i] = buf.mclr[i] & 0x00FFFFFF;
	//if ( buf.mdensity[i] > 0.0014 ) buf.mclr[i] += 0xAA000000;

	// Leap-frog Integration
	float3 vnext = accel*dt + vel;				// v(t+1/2) = v(t-1/2) + a(t) dt	
	//My add  **************************************************************************************************
	//buf.mveleval[i] = (vel + vnext) * 0.5;		// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5			
	//buf.mvel[i] = vnext;
	
	buf.mveleval[i] = (vel + vnext) * 0.5 * phase;		// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5			
	buf.mvel[i] = vnext * phase;
	//**********************************************************************************************************
	buf.mpos[i] += vnext * (dt/ss);						// p(t+1) = p(t) + v(t+1/2) dt		
}


void updateSimParams ( FluidParams* cpufp )
{
	cudaError_t status;
	#ifdef CUDA_42
		// Only for CUDA 4.x or earlier. Depricated in CUDA 5.0+
		// Original worked even if symbol was declared __device__
		status = cudaMemcpyToSymbol ( "simData", cpufp, sizeof(FluidParams) );
	#else
		// CUDA 5.x+. Only works if symbol is declared __constant__
		status = cudaMemcpyToSymbol ( simData, cpufp, sizeof(FluidParams) );
	#endif

	/*app_printf ( "SIM PARAMETERS:\n" );
	app_printf ( "  CPU: %p\n", cpufp );	
	app_printf ( "  GPU: %p\n", &simData );	 */
}


//My add
__device__ int searchIce(int i, float3 ipos, float idens, int cell, bufList buf)
{

	if ( buf.mgridcnt[cell] == 0 ) return 0;	

	int n = 0,j;
	float3 dist;
	float dsq;


	for ( int cndx = buf.mgridoff[ cell ]; cndx < buf.mgridoff[ cell ] + buf.mgridcnt[ cell ]; cndx++ ) {	
		j = buf.mgrid[ cndx ];				
		dist =  ipos - buf.mpos[ j ] ;		// dist in cm   两点间距离
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		if ( dsq < simData.rd2 && dsq > 0) {
			if(buf.mphase[j] < 0.5)
				n++;
		}
	}

	return n;
}

__global__ void updatePhase (bufList buf, int pnumint, float ss)
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if ( i >= pnumint ) return;

	register int newPhase = 1, n = 0;
	register float diff;
	register float3 pos = buf.mpos[i];

	uint gc = buf.mgcell[ i ];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;

	for (int c=0; c < simData.gridAdjCnt; c++) {
		n += searchIce( i, buf.mpos[ i ], buf.mdensity[ i ], gc + simData.gridAdj[c], buf);
	}

	if(n > 15)
	{
		newPhase = 0;
	}

	//diff < 0 说明还没到底，反之说明到底了
	diff = simData.pradius - (pos.y - (simData.pboundmin.y + (pos.x-simData.pboundmin.x)*simData.pground_slope )) * ss;

	if(diff - EPSILON > 0)
	{
		newPhase = 0;
	}

	buf.mphase[i] = newPhase;
}

__device__ float3 contributeTension ( int i, float3 ipos, float idens, int cell, bufList buf)
{
	float dsq, c;	
	float3 dist, tension = make_float3(0,0,0);	
	int j;	

	if ( buf.mgridcnt[cell] == 0 ) return make_float3(0,0,0);	

	tension = make_float3(0,0,0);

	for ( int cndx = buf.mgridoff[ cell ]; cndx < buf.mgridoff[ cell ] + buf.mgridcnt[ cell ]; cndx++ ) {										
		j = buf.mgrid[ cndx ];				
		dist =  ipos - buf.mpos[ j ] ;		// dist in cm   两点间距离
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		if ( dsq < simData.rd2 && dsq > 0) {			
			dsq = dsq * simData.d2;//是实际大小
			c = simData.r2 - dsq ;
			float gradient = c * -6 * simData.poly6kern * simData.r2 * simData.r2;
			float3 dn = simData.psimscale * idens * (buf.mdensity[ j ])*  gradient * dist * (1.0/sqrt(dsq));	
			tension  += dn;
			//buf.msurfaceTension[i] += dn ;
			//buf.msurfaceTension[j] -= 0.5 * dn ;

		}	
	}

	return tension;

}

//__device__ float contributeWeight (int i, float3 ipos, int cell, bufList buf)
//{
//	float total_Weight = 0,dsq,c;
//	int j;
//	float3 dist;
//	for ( int cndx = buf.mgridoff[ cell ]; cndx < buf.mgridoff[ cell ] + buf.mgridcnt[ cell ]; cndx++ ) {										
//		j = buf.mgrid[ cndx ];				
//		dist =  ipos - buf.mpos[ j ];		// dist in cm   两点间距离
//		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
//		if ( dsq < simData.rd2 && dsq > 0) {			
//			c = simData.r2 - dsq * simData.d2 ;
//			total_Weight += c / simData.r2;
//		}
//	}
//	return total_Weight;
//}

__global__ void computeTension(bufList buf, int pnum)
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;

	// Get search cell	
	uint gc = buf.mgcell[ i ];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;

	// 	
	register float3 tension,smoothedtension,n;
	register float total_weight = 0;

	tension = make_float3(0,0,0);
	smoothedtension = make_float3(0,0,0);
	n = make_float3(0,0,0);

	for (int c=0; c < simData.gridAdjCnt; c++) {
		n += contributeTension ( i, buf.mpos[ i ], buf.mdensity[ i ], gc + simData.gridAdj[c], buf);
		//total_weight += contributeWeight ( i, buf.mpos[ i ], gc + simData.gridAdj[c], buf);
	}
	__syncthreads();

	//此时的tension其实是n   loop 2th
	float magnitude = n.x * n.x + n.y * n.y + n.z * n.z;
	magnitude = sqrt(magnitude);

	if(magnitude > 0)
		tension = magnitude * 30 * n;
	buf.msurfaceTension[i] = tension;
	__syncthreads();

	for (int cc=0; cc < simData.gridAdjCnt; cc++) {
		float weight = 0,dsq,c;
		int j;
		int cell = gc + simData.gridAdj[cc];
		float3 dist,ipos;
		for ( int cndx = buf.mgridoff[ cell ]; cndx < buf.mgridoff[ cell ] + buf.mgridcnt[ cell ]; cndx++ ) {										
			j = buf.mgrid[ cndx ];	
			ipos = buf.mpos[i];
			dist =  ipos - buf.mpos[ j ];		// dist in cm   两点间距离
			dsq = dist.x*dist.x + dist.y*dist.y + dist.z*dist.z;
			if ( dsq < simData.rd2 && dsq > 0) {			
				c = simData.r2 - dsq * simData.d2 ;
				weight = c / simData.r2;
				total_weight += weight;
				smoothedtension += weight * buf.msurfaceTension[j];
			}
		}
	}

	if(total_weight > 0)
		buf.mforce[ i ] += smoothedtension/total_weight;
}
