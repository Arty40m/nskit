draw_config = {
    "max_canvas_size":400,
    "nb_colors":{
        "A":"#99627A", 
        "U":"#1F8A70", 
        "T":"#BFDB38", 
        "G":"#E74646", 
        "C":"#068DA9"
    },
    "unknown_nb_color":"#7B8FA1",
    "core_bond_color":"#555", 
    "compl_bond_color":"#526D82", 
    "knot_bond_color":"#B70404"
}


def edit_draw_config(
                max_canvas_size: int = 400,
                A_color: str = "#99627A",
                U_color: str = "#1F8A70",
                T_color: str = "#BFDB38",
                G_color: str = "#E74646",
                C_color: str = "#068DA9",
                unknown_nb_color: str = "#7B8FA1",
                core_bond_color: str = "#555", 
                compl_bond_color: str = "#526D82", 
                knot_bond_color: str = "#B70404"
                     ):
        
        draw_config["max_canvas_size"] = max_canvas_size
        draw_config["nb_colors"]["A"] = A_color
        draw_config["nb_colors"]["U"] = U_color
        draw_config["nb_colors"]["T"] = T_color
        draw_config["nb_colors"]["G"] = G_color
        draw_config["nb_colors"]["C"] = C_color
        draw_config["unknown_nb_color"] = unknown_nb_color
        draw_config["core_bond_color"] = core_bond_color
        draw_config["compl_bond_color"] = compl_bond_color
        draw_config["knot_bond_color"] = knot_bond_color
        