
    # Sub-parser for the itr command
    itr_parser = subparsers.add_parser('itr', help='Itteratively expand a grid around identified minima...')
    itr_parser.add_argument('-path', type=str, default=os.getcwd(), help='Path to read .log files from')
    itr_parser.add_argument('-method', type=int, default=1, help='Energy type to use; 0 = SCF, 1 = counterpoise corrected complexation')
    itr_parser.add_argument('-minima', type=int, default=5, help='Number of minima sites to consider')
    itr_parser.add_argument('-ethr', type=int, default=5, help='Energy cutoff for identifying discrete minima')
    itr_parser.add_argument('-dthr', type=int, default=2, help='Distance cutoff for identifying discrete minima')           
    itr_parser.add_argument('-step', type=int, default=1, help='Iteration step size')
    itr_parser.add_argument('-sweep', type=int, default=1, help='Area to scan around minima (e.g. minima +/- {sweep} Angstroms')
    
	
	
###
# (future) function for iteratively increasing grid size around minima (woah)
###
def handle_itr(args):
# master function
    print(f"Handling iteration with step {args.step}")
    print(f"First, reading data from {args.path}")
    data = process_files(path=args.path, method=args.method)
    dx_values, dy_values, dz_values, dyz_values, e_values, de_values = shape_data(data)
    
    minima_list = find_local_minima(data, num_minima = args.minima, e_threshold = args.ethr , d_threshold = args.dthr)
    count = len(data) # count should be the number of existing files!
    
    for min in minima_list:
        x_new = (min[0] - args.sweep, min[0] + args.sweep)
        y_new = (min[1] - args.sweep, min[1] + args.sweep)
        z_new = (min[2] - args.sweep, min[2] + args.sweep)
        new_grid = geo.gen_grid(x = x_new, y = y_new, z = z_new, res = args.step) # get a new grid around each minima using the user supplied sweep and step values
        outname, count =  create_grid(args.inp, grid = new_grid, x_angle = args.xa, y_angle = args.ya, z_angle = args.za, min_dist = args.min, max_dist = args.max, count = count, filename2 = args.inp2, outname = args.out)
    #
    #gps.make_sge_job(filename = outname, nproc = args.cpu, vmem = args.mem, startjob=1, endjob=count) # make the SGE job
    #
    