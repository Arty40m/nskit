import numpy as np



class CircularGraph:
    
    def calculate_circular_coords(self, n: int):
        N = n + 2
        R = 1/(np.sin(np.pi/N))
        dela = 2*np.pi/N
        rotm = np.array([[np.cos(dela), np.sin(dela)], [-np.sin(dela), np.cos(dela)]], dtype=np.float32)

        vec = np.array([0, -R])@rotm
        vecs = np.zeros((n, 2), dtype=np.float32)
        for i in range(n):
            new_vec = vec@rotm
            nbc = (vec+new_vec)/2
            vecs[i] = nbc
            vec = new_vec
            
        return vecs, R


    def calculate_helix_radiuses(self, nb_coords, R):
        helix_radiuses = []

        for h in self.helixes:
            radiuses = []

            for i, j in h:
                iv = nb_coords[i]
                jv = nb_coords[j]
            
                C = iv + (jv - iv)/2
                clen = np.linalg.norm(C)
                r = (R**2)/clen

                radiuses.append(r)
            helix_radiuses.append(radiuses)

        return helix_radiuses