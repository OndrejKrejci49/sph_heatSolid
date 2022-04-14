#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"
#include "SPH_viscouse_forces.h"
#include "ALESPH_riemann_solver.h"

#include "SPH_pressure.h"


void Compute_Acceleration_BT
(Particle_system &particles)
{

	//Prepare cells to find interacting particles
	unsigned int ncx, ncy;
	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	//Indices of neigbour cells

	#pragma omp parallel for
	for(int c = 0; c < particles.cells.size(); c++)
	{

		idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
		std::vector<int> ac;

		//Load indices of neighbour cells
		zz = c;
		zp = c + 1;
		pp = c + 1 + ncy;
		pz = c + ncy;
		pm = c - 1 + ncy;
		zm = c - 1;
		mm = c - 1 - ncy;
		mz = c - ncy;
		mp = c + 1 - ncy;

		ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "ACCELERATION -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

		//Temp. data, temp. variable
		realvec ac_temp = {0.,0.};
		realvec ac_temp_visco = {0.,0.};
		real p_temp = 0;
		real visco = 0;
		real gamma = 0;
		realvec omegaa_temp = {0., 0.};

		real drs; //dr size
		real drdv; //vect(v) \cdot \nabla W
		realvec dW; //smoothing function gradient
		real W;
		real ndr;
		real rdW;

		real h = particles.data_const.h;
		real m = particles.data_const.m;
		real dp = particles.data_const.dp;
		real c0 = particles.data_const.cs;
		real rho0 = particles.data_const.rho0;

		realvec ar; //position of ative particle
		realvec av; //position of ative particle
		real arho; //actual particle densiy
		real ap; //actual particle density
		real aomega; //ALE

		double *kernel;

		//exprimet for(int i = 0; i < particles.cells[zz].np; i++)
		for(int i = 0; i < particles.cells[zz].cp.size(); i++)
		{

			if(particles.data.part_type[particles.cells[zz].cp[i]] != fluid){continue;}
			//if(particles.data.part_type[particles.cells[zz].cp[i]] == outletf){continue;}
			//if(particles.data.part_type[particles.cells[zz].cp[i]] == outlet){continue;}
			//if(particles.data.part_type[particles.cells[zz].cp[i]] == inlet){continue;}

			/*
			if(particles.data.part_type[particles.cells[zz].cp[i]] == wall || particles.data.part_type[particles.cells[zz].cp[i]] == virt) {continue;}
			*/


			//Load data of actual particle
			ar = particles.data.r[particles.cells[zz].cp[i]];
			av = particles.data.v[particles.cells[zz].cp[i]];
			arho = particles.data.rho[particles.cells[zz].cp[i]];
			ap = particles.data.p[particles.cells[zz].cp[i]];
			aomega = particles.special.omega[particles.cells[zz].cp[i]];


			for(int &cl: ac)
			{

				if( cl < 0 ){continue;}
				//experiment for(int n = 0; n < particles.cells[cl].np; n++)
				for(int n = 0; n < particles.cells[cl].cp.size(); n++)
				{

					if(particles.cells[zz].cp[i] == particles.cells[cl].cp[n]){continue;}
					/*
				if(particles.data.part_type[particles.cells[cl].cp[n]] == outletf){continue;}
				*/
				realvec nr; //position of neighbour particle
				realvec dr; //position diference
				realvec nv; //position of neighbour particle
				realvec dv; //velocity diference
				real nrho; //neighbour density
				real np; //neighbour pressure
				real nomega; //ALE

				//Load data of neighbour particle
				nr = particles.data.r[particles.cells[cl].cp[n]];
				nv = particles.data.v[particles.cells[cl].cp[n]];
				nrho = particles.data.rho[particles.cells[cl].cp[n]];
				np = particles.data.p[particles.cells[cl].cp[n]];
				nomega = particles.special.omega[particles.cells[cl].cp[n]];

				//ALESPH
				real vr, vl; //left and right state
				realvec drn; //unit pair vector
				real rhos; //rho^star, Riemann problem solution
				real vss; //v^star, Riemann problem solution
				realvec vs; //v^star, Riemann problem solution
				real vsdW;
				realvec vz;
				real ps;

				//Position and velocity difference
				dr = ar - nr;
				dv = av - nv;
				drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
				drdv = dr.x*dv.x + dr.y*dv.y;

				/* get kernel values
				double *Wendland_kernel(double r, double h) */
				kernel = Wendland_kernel(drs, h);
				dW = dr * kernel[1]; // <--- check this, if its ok
				W = kernel[0];

				drn = dr*(1/drs)*(-1); //unit pair vector
				vl = av.x * drn.x + av.y * drn.y;
				vr = nv.x * drn.x + nv.y * drn.y;


				rhos = densRiemannLinearized(arho, nrho, vl, vr, 0.5*(arho + nrho), c0);
				vss = velRiemannLinearized(arho, nrho, vl, vr, 0.5*(arho + nrho), c0);
				//vss = velRiemannLinearizedwithPressure(arho, nrho, vl, vr, ap, np, 0.5*(arho + nrho), c0);
				vs = drn * vss + ((av + nv)*0.5 - drn*(vl + vr)*0.5);
				//ps = Compute_Pressure2(rhos,  rho0,  c0);

				//ps = pRiemannLinearized(arho, nrho, vl, vr, ap, np, 0.5*(arho + nrho), c0);
				ps = pRiemannLinearizedWithLimiter(arho, nrho, vl, vr, ap, np, 0.5*(arho + nrho), c0);

				vsdW = (vs.x - dv.x)*dW.x + (vs.y - dv.y)*dW.y;
				vz = (av + nv)*0.5;
				//vz = dv;

				//std::cout << "ACCELERATION >> Star/ rhos: " << rhos <<  " vss: " << vss << " vs: [" << vs.x << "," << vs.y << "]" << " drn: [" << drn.x << "," << drn.y << "]" << " W:" << W << \
				 	"\n 1ID: " << particles.cells[zz].cp[i] << " 2ID: " << particles.cells[cl].cp[n] << " vl: " << vl << " vr:" << vr << " av: [" << av.x << "," << av.y << "]" << " nv: [" << nv.x << "," << nv.y << "]" << std::endl;

				/* pressure term */
				p_temp = ap/pow(arho, 2) + np/pow(nrho,2);
				//p_temp = (ap + np)/(arho*nrho);

				/* real Artificial_Viscosity
				(real h, real drs, real drdv, real rho0, real c0, real alpha) */
				visco = Artificial_Viscosity(h, drs, drdv, particles.data_const.rho0, particles.data_const.cs, particles.data_const.avisc);

				// <--- check this, if its ok

				//ac_temp = ac_temp - dW*(p_temp + visco)*m; //+= operator is not overloaded yet

				if(particles.data.part_type[particles.cells[cl].cp[n]] == wall)
				{

					realvec nb = {0., 0.};
					//nb = particles.special.n[particles.cells[cl].cp[n]]*(-1);
					nb = {0., -1.};
					ndr = dr.x*nb.x + dr.y*nb.y;


					ac_temp.x += p_temp * nb.x * W * dp * nrho;
					ac_temp.y += p_temp * nb.y * W * dp * nrho;

					//omegaa_temp.x -= 2*nomega*W*dp*((rhos*vs.x*(vs.x - vz.x) + ps) * nb.x + (rhos*vs.x*(vs.y - vz.y)) * nb.y);
					//omegaa_temp.y -= 2*nomega*W*dp*((rhos*vs.y*(vs.x - vz.x)) * nb.x + (rhos*vs.y*(vs.y - vz.y) + ps) * nb.y);

				 	ndr = dr.x*nb.x + dr.y*nb.y;

						vl = (-1)*( av.x * nb.x + av.y * nb.y );
						vr = -vl;

						realvec rwf;
						rwf = ar;

						real pR = ap + arho*(-1)*(particles.data_const.graviy.y*rwf.y);

						real vss = velRiemannLinearizedwithPressure(arho, nrho, vl, vr, ap, np, 0.5*(arho + nrho), c0);
						real ps = pRiemannLinearized(arho, nrho, vl, vr, ap, pR, 0.5*(arho + nrho), c0);
						realvec vs = drn * vss + ((av + nv)*0.5 - drn*(vl + vr)*0.5);

						realvec dvs = av - vs;
						real dvsdW = dvs.x*dW.x + dvs.y*dW.y;

				 	real dvn = 0;
				 	//dvn = dv.x*nb.x + dv.y*nb.y; //toto je ok
				 	dvn = dvs.x*nb.x + dvs.y*nb.y;


					omegaa_temp.x -= 2*nomega*W*dp*( ps * nb.x);
					omegaa_temp.y -= 2*nomega*W*dp*( ps * nb.y);


					ac_temp_visco.x = ac_temp_visco.x - dW.x*(visco )*m;
					ac_temp_visco.y = ac_temp_visco.y - dW.y*(visco )*m;

				}
				else
				{

					//rdW = dr.x * dW.x + dr.y * dW.y;

					//ac_temp.x = ac_temp.x - dW.x*(p_temp )*m;
					//ac_temp.y = ac_temp.y - dW.y*(p_temp )*m;

					ac_temp_visco.x = ac_temp_visco.x - dW.x*(visco )*m;
					ac_temp_visco.y = ac_temp_visco.y - dW.y*(visco )*m;

					omegaa_temp.x -= 2*nomega*((rhos*vs.x*(vs.x - vz.x) + ps) * dW.x + (rhos*vs.x*(vs.y - vz.y)) * dW.y);
					omegaa_temp.y -= 2*nomega*((rhos*vs.y*(vs.x - vz.x)) * dW.x + (rhos*vs.y*(vs.y - vz.y) + ps) * dW.y);

				}

				gamma += W*m/nrho;

				/* Temp, step */
				p_temp = 0; visco = 0;

				} // cycle over particles in neighbour cells

			} // cycle over neighbour cells

			if(gamma == 0)
			{

				realvec zero = {0., 0.};
				particles.data.a[particles.cells[zz].cp[i]] = zero - particles.data_const.graviy;
				particles.special.omegaa[particles.cells[zz].cp[i]] = zero - particles.data_const.graviy;

			}
			else
			{

				// particles.data.a[particles.cells[zz].cp[i]] =  ac_temp/gamma +  ac_temp_visco*2*0.001/(gamma * arho) - particles.data_const.graviy;

				//particles.data.a[particles.cells[zz].cp[i]] =  ac_temp/gamma +  ac_temp_visco/gamma - particles.data_const.graviy;
				//particles.data.a[particles.cells[zz].cp[i]] =  ac_temp/gamma   - particles.data_const.graviy;

				particles.data.a[particles.cells[zz].cp[i]] =  ac_temp/gamma +  ac_temp_visco/gamma - particles.data_const.graviy;
				particles.special.omegaa[particles.cells[zz].cp[i]] = omegaa_temp*aomega - particles.data_const.graviy*aomega*arho;

			}

				//std::cout << "ACCELERATION2 >> Write/ aomega: " << particles.special.omega[particles.cells[zz].cp[i]] << std::endl;// " omegaa_temp: " << omegaa_temp <<  std::endl;

				idx index = particles.cells[zz].cp[i];
			//std::cout << "[ACCELERATION2] >> V: [" << particles.data.v[index].x << "," << particles.data.v[index].y  << "] OmegaV: [" << particles.special.omegav[index].x << "," << particles.special.omegav[index].y << "] Omega: " << particles.special.omega[index] << " dOmegaA: [" << particles.special.omegaa[index].x << "," << particles.special.omegaa[index].x << "] dOmega: " << particles.special.domega[index] << " Omega*G: " << particles.data_const.graviy.y*particles.special.omega[index] << std::endl;

			ac_temp = {0., 0.};
			ac_temp_visco = {0., 0.};
			omegaa_temp = {0., 0.};
			gamma = 0;

		} // cycle over particles in active cell

	} // cycle over cells

	#pragma omp barrier

} // function
