#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable


#define GETADRESS(t1,t2,size) ((t1) < (t2)) ? ((size * t1) + t2 - ((t1*(t1+1))/2)) : ((size * t2) + t1 - ((t2*(t2+1))/2))
#define GRID_TAIL_SIZE      64
#define BEAD_NHOOD_SIZE     256

typedef struct Vector3
{
    double x, y, z, w;
} my_Vector3;


typedef struct Bead
{
    my_Vector3    position;
    my_Vector3    image_position;
    my_Vector3    velocity;
    my_Vector3    force;
    my_Vector3    velocity_old;
    my_Vector3    force_old;    
} _Bead;            


typedef struct ParticleInfo
{        
    int          type;
    int          tag;
    int          loc_id;
    double        mass;
    double        mass_inv;
    double        lambda_dt_mass_inv;
    double        half_dt_mass_inv;
    double        half_dt2_mass_inv;
    double        radius;
    double        radius_sq;
    double        skin_radius;
    double        skin_radius_sq;

} _particle_info;



typedef struct GridCellNh
{
    int n[13];
} _gridCellNh;



typedef struct ParticleNh
{
    int n[BEAD_NHOOD_SIZE];
} _particleNh;



typedef struct InteractionNh
{
    double n[BEAD_NHOOD_SIZE];
} _interactionNh;


typedef struct GridCell
{
    int n[GRID_TAIL_SIZE];
} _GridCell;

typedef struct SimulationSettings
{
    int   types_num;
    int   bc_type;
    double wall_contact_stiffness;
    double cutoff;
    double cutoff_sq;
    double box_unit_cell;
    double dt;
    double half_dt;
    double half_dt2;
    double lambda_dt;
    double gamma;
    double sigma;
    double damping;
    int gridCellsX;
    int gridCellsY;
    int gridCellsZ;
    int gridCellsXY;
    double halfBoxSq;
    my_Vector3 gridStep;
    my_Vector3 boxSize;
    my_Vector3 gravity;
} _settings;


typedef struct GridParameters
{
    int gridCellsX;
    int gridCellsY;
    int gridCellsZ;
    int gridCellsXY;
    double halfBoxSq;
    my_Vector3 gridStep;
    my_Vector3 boxSize;
} _parameters;



typedef struct Pair
{
    int type;
    int p1;
    int p2;
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;
    double c6;
} _pair;


typedef struct SpringStructure                                                                    
{     
    int    type;
    int    tag;
    int    status;
    int    p1;
    int    p2;
    double    restLength;
    double    stiffness;
    double    damping;
    double    c1;
    double    c2;                                                                                                                                
} _spring;                                                                        

typedef struct AngleBond
{
    int        type;
    int        tag;
    int        status;
    int        b1;
    int        b2;
    int        b3;
    int        s1;
    int        s2;
    double     angle;
    double     current_angle;
    double      c1;
    double      c2;
    double      c3;
} angle_bond;

inline void atomic_add_global(volatile __global float *addr, float val)
{
       union{
           unsigned int u32;
           float        f32;
       } next, expected, current;
       current.f32    = *addr;
       do{
          expected.f32 = current.f32;
           next.f32     = expected.f32 + val;
           current.u32  = atomic_cmpxchg( (volatile __global unsigned int *)addr, 
                               expected.u32, next.u32);
       } while( current.u32 != expected.u32 );
}


inline void atomic_double_add(volatile __global double *val, double delta)
{
    union 
    {
          double f;
          ulong  i;
          } old;
        
     union
        {
          double f;
          ulong  i;
        } new;
        
    do 
        {
          old.f = *val;
          new.f = old.f + delta;
        } while (atom_cmpxchg ( (volatile __global ulong *)val, old.i, new.i) != old.i);
}

void kernel AdvancePosition(global _Bead* bead_buffer, global _particle_info* bead_info_buffer,  _settings settings, int beads_num)            
{                                                                        
    int bead_id = get_global_id(0);                                                                                                                        
    if(bead_id < beads_num)                                                            
    {    
        _Bead       bead      = bead_buffer[bead_id];
        double       lambda_dt = bead_info_buffer[bead_id].lambda_dt_mass_inv;
        double       half_dt2  = bead_info_buffer[bead_id].half_dt2_mass_inv;
        my_Vector3  increment = bead.velocity_old;

        bead.velocity.x = mad(bead.force.x, lambda_dt, increment.x);                    
        bead.velocity.y = mad(bead.force.y, lambda_dt, increment.y);                
        bead.velocity.z = mad(bead.force.z, lambda_dt, increment.z);

        increment.x *= settings.dt;
        increment.y *= settings.dt;
        increment.z *= settings.dt;

        increment.x = mad(bead.force.x, half_dt2, increment.x);
        increment.y = mad(bead.force.y, half_dt2, increment.y);
        increment.z = mad(bead.force.z, half_dt2, increment.z);

        bead.position.x += increment.x;
        bead.position.y += increment.y;
        bead.position.z += increment.z;

        bead.image_position.x += increment.x;
        bead.image_position.y += increment.y;
        bead.image_position.z += increment.z;
                     
        bead.force_old = bead.force;
        bead.force.x = 0.0;
        bead.force.y = 0.0;
        bead.force.z = 0.0;
                 
 
        bead.force.x = settings.damping * bead.velocity.x;
        bead.force.y = settings.damping * bead.velocity.y;
        bead.force.z = settings.damping * bead.velocity.z;


        bead.force.x += settings.gravity.x;
        bead.force.y += settings.gravity.y;
        bead.force.z += settings.gravity.z;


        if (settings.bc_type == 0)
        {
            if(bead.position.z < 0.0) bead.position.z = settings.boxSize.z;    
            if(bead.position.y < 0.0) bead.position.y = settings.boxSize.y;
            if(bead.position.x < 0.0) bead.position.x = settings.boxSize.x;
            
            if(bead.position.z > settings.boxSize.z) bead.position.z = 0.0;
            if(bead.position.y > settings.boxSize.y) bead.position.y = 0.0;
            if(bead.position.x > settings.boxSize.x) bead.position.x = 0.0;
        }

        if (settings.bc_type == 1)
        {
            if(bead.position.z < 0.0) bead.force.z += bead.position.z * settings.wall_contact_stiffness * -1.0;    
            if(bead.position.y < 0.0) bead.force.y += bead.position.y * settings.wall_contact_stiffness * -1.0; 
            if(bead.position.x < 0.0) bead.force.x += bead.position.x * settings.wall_contact_stiffness * -1.0; 
                
            if(bead.position.z > settings.boxSize.z) bead.force.z += (bead.position.z - settings.boxSize.z) * settings.wall_contact_stiffness * -1.0; 
            if(bead.position.y > settings.boxSize.y) bead.force.y += (bead.position.y - settings.boxSize.y) * settings.wall_contact_stiffness * -1.0; 
            if(bead.position.x > settings.boxSize.x) bead.force.x += (bead.position.x - settings.boxSize.x) * settings.wall_contact_stiffness * -1.0; 
        }
        
        bead_buffer[bead_id] = bead;
    }                                                            
}                                                            

// Resolves angle bond potentials.
void kernel ResolveAngles(global angle_bond* bond, global _spring* springs, global _Bead* bead_buffer, int bond_num)
{
    int bond_id = get_global_id(0);                                                        
    if( bond_id < bond_num)                                                        
    {        
        angle_bond _bond = bond[bond_id];    
	if (_bond.type == 1) 
	{
		if(springs[_bond.s1].status == 0)
		{
			_bond.status = 0;
			bond[bond_id].status = 0;
		}
		if(springs[_bond.s2].status == 0)
		{
			_bond.status = 0;
			bond[bond_id].status = 0;
		}
	}
                                        

	if (_bond.status == 1)
	{
	        my_Vector3 temp = bead_buffer[_bond.b2].image_position;                                                                  
        	my_Vector3 ba = bead_buffer[_bond.b1].image_position; 
	        ba.x -= temp.x;
	        ba.y -= temp.y;
	        ba.z -= temp.z;
        
	        my_Vector3 bc = bead_buffer[_bond.b3].image_position;                                               
	        bc.x -= temp.x;
	        bc.y -= temp.y;
	        bc.z -= temp.z;
        
	        my_Vector3 cb = bc;
	        cb.x*=-1.0;    
	        cb.y*=-1.0;    
	        cb.z*=-1.0;    
        
	        double l_ba = sqrt(ba.x*ba.x + ba.y*ba.y + ba.z*ba.z);
	        double l_bc = sqrt(bc.x*bc.x + bc.y*bc.y + bc.z*bc.z);
                
	        double dotp = (ba.x*bc.x + ba.y*bc.y + ba.z*bc.z) / (l_ba * l_bc);
	        if(dotp < -1.0) dotp = -0.999999;
        	if(dotp >  1.0) dotp =  0.999999;
        
	        double angle = acos(dotp);

	        double magnitude = -2.0 * _bond.c1 * (angle - _bond.angle);
        
	        double fa = magnitude / l_ba;
	        double fc = magnitude / l_bc;
        
        	temp.x = ba.y*bc.z - ba.z*bc.y;
	        temp.y = ba.z*bc.x - ba.x*bc.z;
        	temp.z = ba.x*bc.y - ba.y*bc.x;
               
	        my_Vector3 pa; 
        	my_Vector3 pc; 
        
	        pa.x = (ba.y*temp.z - ba.z*temp.y);
        	pa.y = (ba.z*temp.x - ba.x*temp.z);
	        pa.z = (ba.x*temp.y - ba.y*temp.x);
        
        	double ln = sqrt(pa.x*pa.x + pa.y*pa.y + pa.z*pa.z);
	        if (ln != 0.0) ln = 1.0 / ln;
	        ln *= fa;
	        pa.x *= ln;
	        pa.y *= ln;
	        pa.z *= ln;
                
	        pc.x = (cb.y*temp.z - cb.z*temp.y);
        	pc.y = (cb.z*temp.x - cb.x*temp.z);
	        pc.z = (cb.x*temp.y - cb.y*temp.x);
        
        	ln = sqrt(pc.x*pc.x + pc.y*pc.y + pc.z*pc.z);
	        if (ln != 0.0) ln = 1.0 / ln;
	        ln *= fc;
	        pc.x *= ln;
        	pc.y *= ln;
	        pc.z *= ln;      
        
        	atomic_double_add(&(bead_buffer[_bond.b1].force.x), pa.x);
	        atomic_double_add(&(bead_buffer[_bond.b1].force.y), pa.y);
        	atomic_double_add(&(bead_buffer[_bond.b1].force.z), pa.z);
        
        	atomic_double_add(&(bead_buffer[_bond.b3].force.x), pc.x);
	        atomic_double_add(&(bead_buffer[_bond.b3].force.y), pc.y);
        	atomic_double_add(&(bead_buffer[_bond.b3].force.z), pc.z);
        
	        atomic_double_add(&(bead_buffer[_bond.b2].force.x), -pa.x-pc.x);
		atomic_double_add(&(bead_buffer[_bond.b2].force.y), -pa.y-pc.y);
        	atomic_double_add(&(bead_buffer[_bond.b2].force.z), -pa.z-pc.z);
	}        
    }
}

// Resolves spring/harmonic bond potentials.
void kernel ResolveSprings(global _Bead* bead_buffer, global _spring* springs, int springs_num)
{                                                                    
    int spring_id = get_global_id(0);                                                        
    if( spring_id < springs_num)                                                        
    {        
        _spring     spring  = springs[spring_id]; 

        if (spring.type == 0 || spring.type == 1)
        {        
            if (spring.status == 1)
            {                                                                                                  
                    my_Vector3  ds      = bead_buffer[spring.p1].image_position;                                                            
                    my_Vector3  temp    = bead_buffer[spring.p2].image_position;                                                                                                                    
            
                    ds.x -= temp.x;                                
                    ds.y -= temp.y;                                
                    ds.z -= temp.z;                                

                    double l      = native_sqrt(ds.x*ds.x + ds.y*ds.y + ds.z*ds.z);                        
                    double factor = spring.restLength - l;                                
                    if (l != 0.0) l = 1.0 / l;                                        
                    ds.x *= l;                                                
                    ds.y *= l;                                                
                    ds.z *= l;                                                
                           
                    if (spring.type == 1)
                    {
                    	double stress = -spring.c1 * factor * l;
			if(fabs(stress) > spring.c2)
                        {
                            spring.status = 0;
			    springs[spring_id].status = 0;

			    // bead_buffer[spring.p1].force.x = 0.0;
			    // bead_buffer[spring.p1].force.y = 0.0;
			    // bead_buffer[spring.p1].force.z = 0.0;
			    
			    // bead_buffer[spring.p2].force.x = 0.0;
			    // bead_buffer[spring.p2].force.y = 0.0;
			    // bead_buffer[spring.p2].force.z = 0.0;

			    // bead_buffer[spring.p1].velocity.x = 0.0;
			    // bead_buffer[spring.p1].velocity.y = 0.0;
			    // bead_buffer[spring.p1].velocity.z = 0.0;

			    // bead_buffer[spring.p2].velocity.x = 0.0;
			    // bead_buffer[spring.p2].velocity.y = 0.0;
			    // bead_buffer[spring.p2].velocity.z = 0.0;                               
            

                        }
                    }
	            
		    if (spring.status == 1)
		    {
                    	factor *= spring.stiffness;                                        
	                ds.x*=factor;                                                
        	        ds.y*=factor;                                                
                	ds.z*=factor;        
            
	                atomic_double_add(&(bead_buffer[spring.p1].force.x), ds.x);
        	        atomic_double_add(&(bead_buffer[spring.p1].force.y), ds.y);
                	atomic_double_add(&(bead_buffer[spring.p1].force.z), ds.z);                                
            
            	        ds.x*=-1.0;
                  	ds.y*=-1.0;
	                ds.z*=-1.0;

	                atomic_double_add(&(bead_buffer[spring.p2].force.x), ds.x);
        	        atomic_double_add(&(bead_buffer[spring.p2].force.y), ds.y);
                	atomic_double_add(&(bead_buffer[spring.p2].force.z), ds.z);
		    }
            }
        }
    }                                                    
}                                                    

        
void kernel AdvanceVelocity(global _Bead* bead_buffer, global _particle_info* bead_info_buffer, int beads_num, int buffer_size)
{                                                    
    int bead_id = get_global_id(0);                      
    
    if(bead_id < beads_num)                                        
    {    
            
        double      half_dt      = bead_info_buffer[bead_id].half_dt_mass_inv;
        my_Vector3 velocity_old = bead_buffer[bead_id].velocity_old;
        my_Vector3 force        = bead_buffer[bead_id].force;

        force.x += bead_buffer[bead_id].force_old.x;
        force.y += bead_buffer[bead_id].force_old.y;
        force.z += bead_buffer[bead_id].force_old.z;
        
        velocity_old.x = mad(force.x, half_dt, velocity_old.x);
        velocity_old.y = mad(force.y, half_dt, velocity_old.y);
        velocity_old.z = mad(force.z, half_dt, velocity_old.z);
        bead_buffer[bead_id].velocity_old = velocity_old;   
    }    
}



void kernel ResetBeadEnergy(global _Bead* bead_buffer, int beads_num)
{                                                    
    int bead_id = get_global_id(0);                                        
    if(bead_id < beads_num)                                        
    {    
        my_Vector3 empty_vector;
        empty_vector.x = 0.0;
        empty_vector.y = 0.0;
        empty_vector.z = 0.0;

        bead_buffer[bead_id].force        = empty_vector;
        bead_buffer[bead_id].velocity     = empty_vector;
        bead_buffer[bead_id].force_old    = empty_vector;
        bead_buffer[bead_id].velocity_old = empty_vector;
    }
}



void kernel FillGrid(global _Bead* bead_buffer,  global int* grid_head, global _GridCell* grid_tail, _settings settings, int beads_num)
{                                                    
    int bead_id = get_global_id(0);                                        
    if(bead_id < beads_num)                                        
    {    
        my_Vector3 position = bead_buffer[bead_id].position;
        int cellX = floor((position.x) / settings.gridStep.x);
        int cellY = floor((position.y) / settings.gridStep.y);
        int cellZ = floor((position.z) / settings.gridStep.z);

        if( cellX >= 0 && cellX < settings.gridCellsX &&
            cellY >= 0 && cellY < settings.gridCellsY &&
            cellZ >= 0 && cellZ < settings.gridCellsZ)
        {
            int cell_id = cellX + cellY*settings.gridCellsX + cellZ*settings.gridCellsXY;
            int old = atomic_add(&(grid_head[cell_id]), 1);
            if(old < GRID_TAIL_SIZE) grid_tail[cell_id].n[old] = bead_id;
        }
    }
}

void kernel ResetGrid(global int* grid_head,  int grid_size)
{
    int grid_cell_id = get_global_id(0);
    if (grid_cell_id < grid_size)
    {
        grid_head[grid_cell_id] = 0;
    }
}


void kernel ResetContactList(global int* bead_nh_head, int bead_nh_list_size)
{
    int bead_id = get_global_id(0);
    if (bead_id < bead_nh_list_size)
    {
        bead_nh_head[bead_id] = 0;
    }
}

// Initializes contact list between beads by means of grid search algorithm. 
// Contact is considered as possible when distance between beads is smaller than
// interation cutoff distance. Two consecutive beads from the same chain/fibre cannot
// interact.   

void kernel BuildContactList(global int* grid_head,
                             global _GridCell* grid_tail,
                             global _gridCellNh* grid_cell_nh,
                             global int* bead_nh_head,
                             global _particleNh* bead_nh_tail,
                             global _Bead* bead_buffer,
                             global _particle_info* bead_info,
                             global _particleNh* interaction_id, //nowe
                             __global const _pair* coefficient_table,
                             _settings settings,
                             int grid_size)
{
    int grid_cell_id = get_global_id(0);
    if (grid_cell_id < grid_size)
    {
        bool        is_collidable;
        int         interaction_types_num = settings.types_num;
        double      maxDim  = settings.halfBoxSq;
        double      cutoff_sq = settings.cutoff_sq;
        
        my_Vector3  boxSize = settings.boxSize;
        my_Vector3  rij;
        
        int beads_in_cell = grid_head[grid_cell_id];
        if( beads_in_cell > 0 )
        {
            if(beads_in_cell > GRID_TAIL_SIZE) beads_in_cell = GRID_TAIL_SIZE;
            
            _GridCell    bead_id = grid_tail[grid_cell_id];
            _GridCell    nh_bead_id;
            
            for(int ii = 0; ii < (beads_in_cell-1); ii++)
            {
                int             bead_id_1    = bead_id.n[ii];
                my_Vector3      position_1   = bead_buffer[bead_id_1].position;
                _particle_info  bead_info_1  = bead_info[bead_id_1];
                
                for(int jj = (ii+1); jj < beads_in_cell; jj++)
                {
                    int             bead_id_2    = bead_id.n[jj];
                    _particle_info  bead_info_2  = bead_info[bead_id_2];
                    
                    is_collidable = false;
                    
                    if(bead_info_1.tag != bead_info_2.tag)
                    {
                        is_collidable = true;    
                    }
                    else
                    {
                        if (abs(bead_info_1.loc_id - bead_info_2.loc_id) > 1)
                        {
                            is_collidable = true;
                        }
                    }
                    
                    if(is_collidable == true)
                    { 
                        my_Vector3 position_2 = bead_buffer[bead_id_2].position;
                        
                        rij.x = position_1.x - position_2.x;
                        rij.y = position_1.y - position_2.y;
                        rij.z = position_1.z - position_2.z;
                
                        double length  = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;

                        if (settings.bc_type == 0)
                        {
                             if(length > maxDim)
                             {
                                if(rij.x*rij.x > maxDim)
                                {
                                    if(rij.x > 0.0)
                                        rij.x -= boxSize.x;
                                    else    
                                        rij.x += boxSize.x;
                                }

                                if(rij.y*rij.y > maxDim)
                                {
                                    if(rij.y > 0.0)
                                        rij.y -= boxSize.y;
                                    else    
                                        rij.y += boxSize.y;
                                }

                                if(rij.z*rij.z > maxDim)
                                {
                                    if(rij.z > 0.0)
                                        rij.z -= boxSize.z;
                                    else    
                                        rij.z += boxSize.z;
                                }
                                length = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
                            }
                        }
                        
                        
                        int pair_id = GETADRESS(bead_info_1.type, bead_info_2.type, interaction_types_num);

                        if(coefficient_table[pair_id].type == 1 || coefficient_table[pair_id].type == 3)
                        {
                            if(length < cutoff_sq)
                            {
                                int size = atomic_add(&(bead_nh_head[bead_id_1]), 1);
                                if( size < BEAD_NHOOD_SIZE)
                                {
                                    bead_nh_tail[bead_id_1].n[size] = bead_id_2;
                                    interaction_id[bead_id_1].n[size] = pair_id;
                                }
                            }
                        }
                        
                        if(coefficient_table[pair_id].type == 2)
                        {
                            double min_skin = (bead_info_1.skin_radius > bead_info_2.skin_radius) ? bead_info_2.skin_radius : bead_info_1.skin_radius;
                            double max_radius = (bead_info_1.radius > bead_info_2.radius) ? bead_info_1.radius : bead_info_2.radius;
                            
                            if(length < pow((max_radius + min_skin),2.0))
                            {
                                int size = atomic_add(&(bead_nh_head[bead_id_1]), 1);
                                if( size < BEAD_NHOOD_SIZE)
                                {
                                    bead_nh_tail[bead_id_1].n[size] = bead_id_2;
                                    interaction_id[bead_id_1].n[size] = pair_id;
                                }
                            }
                        }
                        
                        
                    }
                }    
            }
            
            _gridCellNh nh = grid_cell_nh[grid_cell_id];
            for(int ii = 0; ii < 13; ii++)    
            {
                int nh_cell_id = nh.n[ii];
                if( nh_cell_id != -1 )
                {
                    int nh_beads_in_cell = grid_head[nh_cell_id];
                    if( nh_beads_in_cell > 0 )
                    {
                        nh_bead_id = grid_tail[nh_cell_id];
                        
                        if( nh_beads_in_cell > GRID_TAIL_SIZE ) nh_beads_in_cell = GRID_TAIL_SIZE;
                        for(int jj = 0; jj < nh_beads_in_cell; jj++)
                        {
                            int             bead_id_2    = nh_bead_id.n[jj];
                            _particle_info  bead_info_2  = bead_info[bead_id_2];
                            my_Vector3      position_2   = bead_buffer[bead_id_2].position;
                            
                            for(int kk = 0; kk < (beads_in_cell); kk++)
                            {
                                int             bead_id_1    = bead_id.n[kk];
                                _particle_info  bead_info_1  = bead_info[bead_id_1];
                                
                                
                                is_collidable = false;
                    
                                if(bead_info_1.tag != bead_info_2.tag)
                                {
                                    is_collidable = true;    
                                }
                                else
                                {
                                    if (abs(bead_info_1.loc_id - bead_info_2.loc_id) > 1)
                                    {
                                        is_collidable = true;
                                    }
                                }
                                
                                if(is_collidable == true)
                                {
                                    my_Vector3 position_1 = bead_buffer[bead_id_1].position;
                                    rij.x = position_1.x - position_2.x;
                                    rij.y = position_1.y - position_2.y;
                                    rij.z = position_1.z - position_2.z;
                                    
                                    double length  = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
                                    
                                    if (settings.bc_type == 0)
                                    {
                                         if(length > maxDim)
                                         {
                                            if(rij.x*rij.x > maxDim)
                                            {
                                                if(rij.x > 0.0)
                                                    rij.x -= boxSize.x;
                                                else    
                                                    rij.x += boxSize.x;
                                            }

                                            if(rij.y*rij.y > maxDim)
                                            {
                                                if(rij.y > 0.0)
                                                    rij.y -= boxSize.y;
                                                else    
                                                    rij.y += boxSize.y;
                                            }

                                            if(rij.z*rij.z > maxDim)
                                            {
                                                if(rij.z > 0.0)
                                                    rij.z -= boxSize.z;
                                                else    
                                                    rij.z += boxSize.z;
                                            }
                                                
                                            length = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
                                        }
                                    }
                                    
                                    int pair_id = GETADRESS(bead_info_1.type, bead_info_2.type, interaction_types_num);
                                    
                                    if(coefficient_table[pair_id].type == 1 || coefficient_table[pair_id].type == 3)
                                    {
                                        if(length < cutoff_sq)
                                        {
                                            int size = atomic_add(&(bead_nh_head[bead_id_1]), 1);
                                            if( size < BEAD_NHOOD_SIZE)
                                            {
                                               bead_nh_tail[bead_id_1].n[size] = bead_id_2;
                                               interaction_id[bead_id_1].n[size] = pair_id;
                                            }
                                        }
                                    }
                                    
                                    
                                    if(coefficient_table[pair_id].type == 2)
                                    {
                                        double min_skin = (bead_info_1.skin_radius > bead_info_2.skin_radius) ? bead_info_2.skin_radius : bead_info_1.skin_radius;
                                        double max_radius = (bead_info_1.radius > bead_info_2.radius) ? bead_info_1.radius : bead_info_2.radius;
                                        
                                        if(length < pow((max_radius + min_skin),2.0))
                                        {
                                            int size = atomic_add(&(bead_nh_head[bead_id_1]), 1);
                                            if( size < BEAD_NHOOD_SIZE)
                                            {
                                               bead_nh_tail[bead_id_1].n[size] = bead_id_2;
                                               interaction_id[bead_id_1].n[size] = pair_id;
                                            }
                                        }
                                    }
                                   
                                   
                                }
                            }
                        }
                    }
                }
            }
      
      
           
        }
    }
}


void kernel ResolveContacts(global int* bead_nh_head,
                            global _particleNh* bead_nh_tail,
                            global _Bead* bead_buffer,
                            global _particle_info* bead_info,
                            global _particleNh* interaction_id,
                            __global const _pair* coefficient_table,
                            _settings settings,
                            int beads_num)
{
    int bead_id_1 = get_global_id(0);
    
    if (bead_id_1 < beads_num)
    {
        int number_of_neigbours = bead_nh_head[bead_id_1];
        
        if(number_of_neigbours > 0)
        {
            if(number_of_neigbours > BEAD_NHOOD_SIZE) number_of_neigbours = BEAD_NHOOD_SIZE;
            
            bool   has_contact = false; 
            double cutoff_sq = settings.cutoff_sq;
            int    pair_id   = -1;
            double maxDim    = settings.halfBoxSq;
            double radius_1  = bead_info[bead_id_1].radius;
            
            my_Vector3 AccForce;
            AccForce.x = 0.0;
            AccForce.y = 0.0;
            AccForce.z = 0.0;
            
            my_Vector3 Reaction;
            Reaction.x = 0.0;
            Reaction.y = 0.0;
            Reaction.z = 0.0;

            my_Vector3 boxSize    = settings.boxSize;
            my_Vector3 position_1 = bead_buffer[bead_id_1].position;

            my_Vector3 position_2;
            my_Vector3 rij;

            for(int ii = 0; ii < number_of_neigbours; ii++)
            {   
                int bead_id_2 = bead_nh_tail[bead_id_1].n[ii];
                position_2 = bead_buffer[bead_id_2].position;

                rij.x = position_1.x - position_2.x;
                rij.y = position_1.y - position_2.y;
                rij.z = position_1.z - position_2.z;
    
                double length  = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;

                if(length > maxDim)
                {
                    if(rij.x*rij.x > maxDim)
                    {
                        if(rij.x > 0.0)
                            rij.x -= boxSize.x;
                        else    
                            rij.x += boxSize.x;
                     }

                    if(rij.y*rij.y > maxDim)
                    {
                        if(rij.y > 0.0)
                            rij.y -= boxSize.y;
                        else    
                            rij.y += boxSize.y;
                    }

                    if(rij.z*rij.z > maxDim)
                    {
                        if(rij.z > 0.0)
                            rij.z -= boxSize.z;
                        else    
                            rij.z += boxSize.z;
                    }
                    length = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
                }
       
                pair_id = interaction_id[bead_id_1].n[ii];

		//LJ interaction potential
                if(coefficient_table[pair_id].type == 1)
                {
                    if(length < cutoff_sq) 
                    { 
                        double r2      = length;
                        length = native_sqrt(length);
                            
                        if (length < coefficient_table[pair_id].c3)
                        {
                            double sigma2  = coefficient_table[pair_id].c1; 
                            double epsilon = coefficient_table[pair_id].c2;
                            sigma2 *= sigma2;

                            double fr2 = sigma2 / r2;
                            double fr6 = fr2 * fr2 * fr2; 
                            double fpr = 48.0 * epsilon * fr6 * (fr6 - 0.5) / r2;
                            
                            fpr -= coefficient_table[pair_id].c4;
                              
                            if(length != 0.0)
                                length = 1.0/length;
                            else
                                length = 0.0;

                            rij.x *= length;
                            rij.y *= length;
                            rij.z *= length;

                            Reaction.x = fpr * rij.x;
                            Reaction.y = fpr * rij.y;
                            Reaction.z = fpr * rij.z;
                              
                           
                            AccForce.x += Reaction.x;
                            AccForce.y += Reaction.y;
                            AccForce.z += Reaction.z;
                            
                            atomic_double_add(&(bead_buffer[bead_id_2].force.x), -Reaction.x);
                            atomic_double_add(&(bead_buffer[bead_id_2].force.y), -Reaction.y);
                            atomic_double_add(&(bead_buffer[bead_id_2].force.z), -Reaction.z);
                            
                            has_contact = true;
                        }
                    }
                }

		//penalty contact model                
                if(coefficient_table[pair_id].type == 2)
                {
                    length = native_sqrt(length);
                        
                    double radius_2 = bead_info[bead_id_2].radius;
                    double K = coefficient_table[pair_id].c1;
                    double penetration = (radius_1 + radius_2) - length;
                        
                    if (penetration > 0.0)
                    {
                        if(length != 0.0)
                            length = 1.0/length;
                        else
                            length = 0.0;

                        rij.x *= length;
                        rij.y *= length;
                        rij.z *= length;
                            
                        penetration *= K;
                            
                        Reaction.x = penetration * rij.x;
                        Reaction.y = penetration * rij.y;
                        Reaction.z = penetration * rij.z;
                        
                        AccForce.x += Reaction.x;
                        AccForce.y += Reaction.y;
                        AccForce.z += Reaction.z;

                        atomic_double_add(&(bead_buffer[bead_id_2].force.x), -Reaction.x);
                        atomic_double_add(&(bead_buffer[bead_id_2].force.y), -Reaction.y);
                        atomic_double_add(&(bead_buffer[bead_id_2].force.z), -Reaction.z);
                        
                        has_contact = true;
                    }        
                }
                    
                //two-step harmonic potential interaction
                if(coefficient_table[pair_id].type == 3)
                {
                    if (length < coefficient_table[pair_id].c4)
                    {
                        length = native_sqrt(length);
                        double K = coefficient_table[pair_id].c1;
                        double l0 = coefficient_table[pair_id].c2;
                        double factor = (l0 - length);
			
			
			if(factor < coefficient_table[pair_id].c3)
				factor *= K;
			else
				factor *= K;
			        factor += (coefficient_table[pair_id].c6 - length) * coefficient_table[pair_id].c5;

                        if(length != 0.0)
                            length = 1.0/length;
                        else
                            length = 0.0;

                        rij.x *= length;
                        rij.y *= length;
                        rij.z *= length;
                                 
                        Reaction.x = factor * rij.x;
                        Reaction.y = factor * rij.y;
                        Reaction.z = factor * rij.z;
                        
                        AccForce.x += Reaction.x;
                        AccForce.y += Reaction.y;
                        AccForce.z += Reaction.z;

                        atomic_double_add(&(bead_buffer[bead_id_2].force.x), -Reaction.x);
                        atomic_double_add(&(bead_buffer[bead_id_2].force.y), -Reaction.y);
                        atomic_double_add(&(bead_buffer[bead_id_2].force.z), -Reaction.z);
                        
                        has_contact = true;
                    }    
                }    
            }
            
            if(has_contact == true)
            {
                atomic_double_add(&(bead_buffer[bead_id_1].force.x), AccForce.x);
                atomic_double_add(&(bead_buffer[bead_id_1].force.y), AccForce.y);
                atomic_double_add(&(bead_buffer[bead_id_1].force.z), AccForce.z);
            }

        }
    }
}


    