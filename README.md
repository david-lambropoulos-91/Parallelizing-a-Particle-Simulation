# Parallelizing-a-Particle-Simulation

**Serial Algorithm**


    for step: 0->NSTEPS
        
        // Walk through bins apply force on the particles from adjacent bins and to adjacent bins
        for i: 0->#bins
            for j: 0->#bins
                Grab bin at (i,j)
                Set x,y acceleration of particles in bin (i,j) to 0
                for dx: -1->1
                    for dy: -1->1
                        if neighbor is within grid (vertically+horizontally)
                            Grab neighbor bin at (i+dx, j+dy)
                            for k: 0-># particles in bin (i,j)
                                for l: 0-># particles in bin (i,j)
                                    apply_force(bin(i,j)[k], bin(i+dx,j+dy)[l])
        
        // Move respective particles based on forces attacted upon them
        for i: 0->#bins
            for j: 0->#bins
                Grab bin at (i,j)
                for k: 0-># particles in bin (i,j)
                    move(bin(i,j)[k])
                    x = bin(i,j)[k].x/binSize
                    y = bin(i,j)[k].y/binSize
                    if x and y are in bin (i,j)
                        k++
                    else
                        store particle from bin(i,j)[k] to temporary storage
                        remove particle from bin(i,j)[k]
                resize bin(i,j) to size k
        
        // Adjust particles that are no longer in their original bin
        for i: 0->#particles changed bins
            x = particle[i].x/binSize
            y = particle[i].y/binSize
            add particle from temporary storage to bin(x,y)
        empty temporary storage
