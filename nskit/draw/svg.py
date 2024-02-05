import numpy as np
import math

from .config import draw_config



class DrawSVG:    
    
    def make_svg(self, nb_coords, helix_radiuses, R):

        canvas = min(draw_config["max_canvas_size"], 150 + (250/80)*len(self))
        svg = [('<svg version = "1.1" '
                'baseProfile="full" '
                'xmlns = "http://www.w3.org/2000/svg" ' 
                'xmlns:xlink = "http://www.w3.org/1999/xlink" '
                'xmlns:ev = "http://www.w3.org/2001/xml-events" '
                f'height = "{canvas}px"  width = "{canvas}px">')
                ]
        
        r = R+1
        nb_diam = canvas/(2*r)
        nb_coords = nb_diam*(nb_coords*np.array([1,-1]) + r)
        helix_radiuses = [[r*nb_diam for r in h] for h in helix_radiuses]
        
        stroke = max(1.5, 1.5*(canvas/100)*(15/(len(self)+3)))
        radius = (nb_diam - stroke/2)*(25/(stroke+30) + 0.1)

        # nucleic bases
        nucleotides_lines = self._draw_nucleotides(nb_coords, stroke, radius)
        svg.extend(nucleotides_lines)

        # core bonds
        core_bonds_lines = self._draw_core_bonds(nb_coords, stroke, radius)
        svg.extend(core_bonds_lines)

        # complementary bonds
        compl_bonds_lines = self._draw_compl_bonds(nb_coords, helix_radiuses, stroke, radius)
        svg.extend(compl_bonds_lines)

        svg.append("</svg>")
        svg = '\n'.join(svg)

        return svg
    

    def _draw_nucleotides(self, nb_coords, stroke, radius):
        nb_fontsize = max(1, 2.5*(radius - stroke/2)/math.sqrt(2))
        text_shift = 0.35*nb_fontsize
        text_coords = nb_coords + np.array([-text_shift, text_shift])
        
        lines = [f'<g stroke-width="{stroke:.2f}" fill="none">']
        text_lines = [f'<g font-weight="bold" font-size="{nb_fontsize:.2f}px">']

        for i, nb in enumerate(self.seq):
            color = draw_config["nb_colors"].get(nb, draw_config["unknown_nb_color"])
            x, y = nb_coords[i]
            tx, ty = text_coords[i]
            lines.append(f'    <circle cx="{x:.1f}px" cy="{y:.1f}px" r="{radius:.1f}px" stroke="{color}"/>')
            text_lines.append(f'    <text x="{tx:.1f}" y="{ty:.1f}">{nb}</text>')
            
        lines.append('</g>')
        text_lines.append('</g>')
        lines.extend(text_lines)

        return lines
    

    def _draw_core_bonds(self, nb_coords, stroke, radius):
        lines = [f'<g stroke="{draw_config["core_bond_color"]}" stroke-width="{stroke:.2f}">']
        shift = radius + stroke/2
        
        for j in range(1, nb_coords.shape[0]):
            iv = nb_coords[j-1]
            jv = nb_coords[j]
            ijv = (ijv:=jv-iv)/np.linalg.norm(ijv)
            
            d1 = iv + ijv*shift
            d2 = jv - ijv*shift
            lines.append(f'    <line x1="{d1[0]:.1f}" y1="{d1[1]:.1f}" x2="{d2[0]:.1f}" y2="{d2[1]:.1f}"/>')
        
        lines.append('</g>')
        return lines
    

    def _draw_compl_bonds(self, nb_coords, helix_coords, stroke, radius):
        lines = [f'<g stroke-width="{stroke:.2f}" fill="none" >']
        shift = radius + stroke/2
        knots = self.knots
        
        for n, (h, hr) in enumerate(zip(self.helixes, helix_coords)):
            color = draw_config['knot_bond_color'] if n in knots else draw_config['compl_bond_color']

            for i, (nb1, nb2) in enumerate(h):
                if abs(nb2-nb1+1)>((len(self)+2)//2):
                    nb1, nb2 = nb2, nb1

                iv = nb_coords[nb1]
                r = hr[i]
                jv = nb_coords[nb2]

                ijv = (ijv:=jv-iv)/np.linalg.norm(ijv)
                d1 = iv + ijv*shift
                d2 = jv - ijv*shift
                arc = f'    <path  d="M {d1[0]:.1f} {d1[1]:.1f} A {r:.1f} {r:.1f} 0 0 1 {d2[0]:.1f} {d2[1]:.1f}" stroke="{color}"/>'
                lines.append(arc)
        
        lines.append('</g>')
        return lines