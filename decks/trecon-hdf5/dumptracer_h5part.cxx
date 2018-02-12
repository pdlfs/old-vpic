{
	char fname[256];

	h5part_int64_t ierr;
	int np_local,ngrid[3];
  float location[6];
	species_t *sp; 

	H5PartFile * h5pf;

	h5part_float32_t *Pf ;
	h5part_int32_t *Pi ;

	// get the total number of particles. in this example, output only electrons
	//sp = find_species ("electron_tracer");
  sp = global->tracers_list;

  while(sp){
	np_local = sp->np; // number of particles on this rank

	Pf = (h5part_float32_t *) sp->p;
	Pi = (h5part_int32_t *) sp->p;

	// open H5part file in "particle/T.<step>/" subdirectory
	// filename: eparticle.h5p
	//sprintf (fname, "tracer/T.%d", step);
	//dump_mkdir(fname);
	sprintf (fname, "tracer/T.%d/%s.h5p", step, sp->name);

	double el1 = mp_elapsed(grid->mp);
	//h5pf = H5PartOpenFileParallel (fname, H5PART_WRITE | H5PART_FS_LUSTRE, MPI_COMM_WORLD);
	h5pf = H5PartOpenFileParallel (fname, H5PART_WRITE, MPI_COMM_WORLD);

	ierr = H5PartSetStep (h5pf, step);  
	ierr = H5PartSetNumParticlesStrided (h5pf, np_local, 8);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in setting the # of tracers; Error code: " << ierr);

	//el1 = mp_elapsed(grid->mp) - el1;
	//sim_log("Time in opening and setting up H5part file for tracers: "<< el1 << " s"); 

	//double el2 = mp_elapsed(grid->mp);

	ierr = H5PartWriteDataFloat32 (h5pf, "dX", Pf);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dX in tracer; Error code: " << ierr);
	//sim_log("Finished writing dX");

	ierr = H5PartWriteDataFloat32 (h5pf, "dY", Pf+1);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dY in tracer; Error code: " << ierr);
	//sim_log("Finished writing dY");

	ierr = H5PartWriteDataFloat32 (h5pf, "dZ", Pf+2);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dZ in tracer; Error code: " << ierr);
	//sim_log("Finished writing dZ");

	ierr = H5PartWriteDataInt32   (h5pf, "i",  Pi+3);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing i in tracer; Error code: " << ierr);
	//sim_log("Finished writing i");

	ierr = H5PartWriteDataFloat32 (h5pf, "Ux", Pf+4);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing Ux in tracer; Error code: " << ierr);
	//sim_log("Finished writing Ux");

	ierr = H5PartWriteDataFloat32 (h5pf, "Uy", Pf+5);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing Uy in tracer; Error code: " << ierr);
	//sim_log("Finished writing Uy");

	ierr = H5PartWriteDataFloat32 (h5pf, "Uz", Pf+6);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing Uz in tracer; Error code: " << ierr);
	//sim_log("Finished writing Uz");

	//ierr = H5PartWriteDataFloat32 (h5pf, "q", Pf+7);
	ierr = H5PartWriteDataInt32 (h5pf, "q", Pi+7);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing q in tracer; Error code: " << ierr);
	//sim_log("Finished writing q");

	//el2 = mp_elapsed(grid->mp) - el2;
	el1 = mp_elapsed(grid->mp) - el1;
	sim_log("Time in writing H5part file for tracer: "<< el1 << " s"); 

	double el3 = mp_elapsed(grid->mp);
	ierr = H5PartCloseFile (h5pf);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in closing " << fname << ". Error code: " << ierr);

	el3 = mp_elapsed(grid->mp) - el3;
	sim_log("Time in closing H5part file for tracer: "<< el3 << " s"); 

	// Write metadata if step == 0
	char meta_fname[256];

	sprintf (meta_fname, "tracer/T.%d/grid_metadata_%s.h5p", step, sp->name);

	double el4 = mp_elapsed(grid->mp);
	H5PartFile * h5pmf;

	h5pmf = H5PartOpenFileParallel (meta_fname, H5PART_WRITE, MPI_COMM_WORLD);

	h5part_float32_t *Px0, *Py0, *Pz0, *Pdx, *Pdy, *Pdz ;
	h5part_int32_t *Pnp_local, *Pnx, *Pny, *Pnz ;

	ngrid[0] =  grid->nx;
	ngrid[1] =  grid->ny;
	ngrid[2] =  grid->nz;
	location[0] =  grid->x0;
	location[1] =  grid->y0;
	location[2] =  grid->z0;
	location[3] =  grid->dx;
	location[4] =  grid->dy;
	location[5] =  grid->dz;

	ierr = H5PartSetStep (h5pmf, step);
	//ierr = H5PartSetNumParticlesStrided (h5pmf, 1, 10);
	ierr = H5PartSetNumParticles (h5pmf, 1);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in setting the # of grid metadata entries; Error code: " << ierr);

  ierr = H5PartWriteDataInt32(h5pmf, "np_local", (h5part_int32_t *) &np_local);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing np_local; Error code: " << ierr);
  ierr = H5PartWriteDataInt32(h5pmf, "ngrid", (h5part_int32_t *) &ngrid);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing np_local; Error code: " << ierr);
  ierr = H5PartWriteDataFloat32(h5pmf, "location", (h5part_float32_t *) &location);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing np_local; Error code: " << ierr);


	//sim_log("Finished writing np_local");
  //ierr = H5PartWriteDataInt32(h5pmf, "nx", (h5part_int32_t *) &grid->nx);
  //if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing nx; Error code: " << ierr);
  //sim_log("Finished writing nx");
  //
  //ierr = H5PartWriteDataInt32(h5pmf, "ny", (h5part_int32_t *) &grid->ny);
  //if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing ny; Error code: " << ierr);
  //sim_log("Finished writing ny");
 
  //ierr = H5PartWriteDataInt32(h5pmf, "nz", (h5part_int32_t* ) &grid->nz);
  //if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing nz; Error code: " << ierr);
  //sim_log("Finished writing nz");
  
  //ierr = H5PartWriteDataFloat32(h5pmf, "x0", (h5part_float32_t* ) &grid->x0);
  //if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing x0; Error code: " << ierr);
  //sim_log("Finished writing x0");
  //
  //ierr = H5PartWriteDataFloat32(h5pmf, "y0", (h5part_float32_t* ) &grid->y0);
  //if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing y0; Error code: " << ierr);
  //sim_log("Finished writing y0");
 
 //ierr = H5PartWriteDataFloat32(h5pmf, "z0", (h5part_float32_t* ) &grid->z0);
  //if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing z0; Error code: " << ierr);
  //sim_log("Finished writing z0");
  //
  //ierr = H5PartWriteDataFloat32(h5pmf, "dx", (h5part_float32_t* ) &grid->dx);
  //if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dx; Error code: " << ierr);
  //sim_log("Finished writing dx");
  //
  //ierr = H5PartWriteDataFloat32(h5pmf, "dy", (h5part_float32_t* ) &grid->dy);
  //if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dy; Error code: " << ierr);
  //sim_log("Finished writing dy");
  //
  //ierr = H5PartWriteDataFloat32(h5pmf, "dz", (h5part_float32_t* ) &grid->dz);
  //if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dz; Error code: " << ierr);
  //sim_log("Finished writing dz");
  //



	ierr = H5PartCloseFile (h5pmf);
	if (ierr != H5PART_SUCCESS) sim_log ("Error occured in closing " << meta_fname << ". Error code: " << ierr);
	el4 = mp_elapsed(grid->mp) - el4;
	sim_log("Time in writing metadata: "<< el4 << " s"); 

  sp = sp->next;

  }
}
