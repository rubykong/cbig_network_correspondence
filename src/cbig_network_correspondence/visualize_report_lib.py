import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import colorcet as cc
from matplotlib.lines import Line2D
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from plottable import Table, ColumnDefinition
from plottable.formatters import decimal_to_percent
from plottable.plots import bar, percentile_bars, percentile_stars, progress_donut
from plottable.cmap import normed_cmap
from plottable.plots import circled_image # image
from pandas.plotting import table
import os
from . import visualize_overlap_lib as vis

def dice_percent(val: float) -> str:
    """Formats Numbers to a string, replacing
        0 with "–"
        1 with "✓"
        values < 0.01 with "<1%" and
        values > 0.99 with ">99%"

    Args:
        val (float): numeric value to format

    Returns:
        str: formatted numeric value as string
    """
    if val == 0:
        return "–"
    elif val < 0.01:
        return "<1%"
    elif val > 0.99:
        return ">99%"
    else:
        return f"{str(round(val * 100))}%"

def pvalue_format(val: float) -> str:
    """Formats Numbers to a string, replacing
        0 with "–"
        1 with "✓"
        values < 0.01 with "<1%" and
        values > 0.99 with ">99%"

    Args:
        val (float): numeric value to format

    Returns:
        str: formatted numeric value as string
    """
    if val == 0:
        return "–"
    elif val < 0.05:
        # add * at the end of the string
        return f"{str(round(val,4))}*"
    else:
        # return a float with 4 digits after decimal point
        return f"{str(round(val,4))}"

def dice_format(val: float) -> str:
    """Formats Numbers to a string, replacing
        0 with "–"
        1 with "✓"
        values < 0.01 with "<1%" and
        values > 0.99 with ">99%"

    Args:
        val (float): numeric value to format

    Returns:
        str: formatted numeric value as string
    """
    if val == 0:
        return "–"
    else:
        # return a float with 4 digits after decimal point
        return f"{str(round(val,4))}"



def get_label_rotation(angle, offset):
    # Rotation must be specified in degrees :(
    rotation = np.rad2deg(angle + offset)
    if angle <= np.pi:
        alignment = "right"
        rotation = rotation + 180
    else: 
        alignment = "left"
    return rotation, alignment

def get_group_rotation(angle):
    # Rotation must be specified in degrees :(
    rotation = np.rad2deg(angle)
    if angle <= np.pi:
        alignment = "left"
        rotation = rotation + 270
    else: 
        alignment = "right"
        rotation = rotation + 90

    return rotation, alignment

def add_labels(angles, values, labels, offset, ax, fontsizes):
    
    # This is the space between the end of the bar and the label
    padding = 6
    
    # Iterate over angles, values, and labels, to add all of them.
    for angle, label, fontsize in zip(angles, labels,fontsizes):
        angle = angle
        
        # Obtain text rotation and alignment
        rotation, alignment = get_label_rotation(angle, offset)

        # And finally add the text
        ax.text(
            x=angle, 
            y=values + padding, 
            s=label, 
            ha=alignment, 
            va="center",
            rotation=rotation, 
            rotation_mode="anchor",
            size=fontsize
        )

def add_labels_equalsize(angles, values, labels, offset, ax, fontsize):
    
    # This is the space between the end of the bar and the label
    padding = 6
    
    # Iterate over angles, values, and labels, to add all of them.
    for angle, label in zip(angles, labels):
        angle = angle
        
        # Obtain text rotation and alignment
        rotation, alignment = get_label_rotation(angle, offset)

        # And finally add the text
        ax.text(
            x=angle, 
            y=values + padding, 
            s=label, 
            ha=alignment, 
            va="center",
            rotation=rotation, 
            rotation_mode="anchor",
            size=fontsize,
            alpha=0.8
        ) 

def generate_random_color():
    # Generate a random RGB color
    color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
    return color

def scale_to_interval(x, low=1, high=50):
    MAX = np.max(x)
    MIN = np.min(x)
    sc = ((x - MIN) / (MAX - MIN)) * (high - low) + low
    sc = sc.astype(np.float64)
    return sc

def scale_fontsize_uni(x, low=14, high=18):
    MAX = np.max(x)
    MIN = np.min(x)
    sc = ((x - MIN) / (MAX - MIN)) * (high - low) + low
    sc = sc.astype(np.float64)
    return sc

def scale_fontsize(x, groups, low=12, high=20):
    #extract min and max for each group
    MAX = []
    MIN = []
    sc = np.zeros(x.shape)
    for group in groups:
        MAX = np.max(x[groups==group])
        MIN = np.min(x[groups==group])
        sc[groups==group] = ((x[groups==group] - MIN) / (MAX - MIN)) * (high - low) + low
        
    sc = sc.astype(np.float64)
    return sc

def circular_barchart(df,group_names,group_names_full,p_threshold=0.05):
    # if group_names is a list, change it to numpy array
    if type(group_names) == list:
        group_names = np.array(group_names)
    # if group_names_full is a list, change it to numpy array
    if type(group_names_full) == list:
        group_names_full = np.array(group_names_full)
    # All this part is like the code above
    for group in group_names:
        curr_VALUES = df[df["group"]==group]["dice"].values * 100
        curr_LABELS = df[df["group"]==group]["name"].values
        curr_GROUP = df[df["group"]==group]["group"].values
        curr_pvalue = df[df["group"]==group]["p_value"].values
    
        # concatente all values and labels
        if group == group_names[0]:
            VALUES = curr_VALUES
            LABELS = curr_LABELS
            GROUP = curr_GROUP
            PVALUE = curr_pvalue
        else:
            VALUES = np.concatenate((VALUES,curr_VALUES))
            LABELS = np.concatenate((LABELS,curr_LABELS))
            GROUP = np.concatenate((GROUP,curr_GROUP))
            PVALUE = np.concatenate((PVALUE,curr_pvalue))
            
    # if p_value is higher than p_threshold, set Label to ""
    LABELS[PVALUE>p_threshold] = ""

    PAD = 3
    ANGLES_N = len(VALUES) + PAD * len(np.unique(GROUP))
    ANGLES = np.linspace(0, 2 * np.pi, num=ANGLES_N, endpoint=False)
    #WIDTH = (2 * np.pi) / len(ANGLES)

    GROUPS_SIZE = [len(df[df["group"]==i]) for i in group_names]

    PLUS = 50
    OFFSET = np.pi / 2
    offset = 0
    IDXS = []
    for size in GROUPS_SIZE:
        IDXS += list(range(offset + PAD, offset + size + PAD))
        offset += size + PAD

    fig, ax = plt.subplots(figsize=(12, 13), subplot_kw={"projection": "polar"})
    fig.patch.set_facecolor('white')

    ax.set_theta_offset(OFFSET)
    ax.set_ylim(0, 100 + PLUS)
    ax.set_frame_on(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    
    palette = sns.color_palette(cc.glasbey, n_colors=len(group_names))

    COLORS = [palette[i] for i, size in enumerate(GROUPS_SIZE) for _ in range(size)]

    ax.vlines(
        ANGLES[IDXS], PLUS,VALUES +  PLUS, color=COLORS, linewidth=1
        #edgecolor="white", linewidth=2
    )
    ax.scatter(ANGLES[IDXS], VALUES +  PLUS , s=scale_to_interval(VALUES),color=COLORS)
    
    add_labels(ANGLES[IDXS], np.max(VALUES) + PLUS, LABELS, OFFSET, ax,scale_fontsize_uni(VALUES))

    # Extra customization below here --------------------

    # This iterates over the sizes of the groups adding reference
    # lines and annotations.

    offset1 = 0 
    for group, size in zip(group_names_full, GROUPS_SIZE):
        #get idx of group, group_names_full is a numpy array
        idx = np.where(group_names_full==group)[0][0] + 1
        # Add line below bars
        x1 = np.linspace(ANGLES[offset1 + PAD], ANGLES[offset1 + size + PAD - 1], num=50)
        ax.plot(x1, [PLUS-10] * 50, color="#333333",linewidth=2)
        x1_rotation, x1_alignment = get_group_rotation(np.mean(x1))
        # Add text to indicate group
        ax.text(
            np.mean(x1), PLUS-20, idx, color="#333333", fontsize=16, 
            fontweight="bold", 
            #ha=x1_alignment, va="center",rotation=x1_rotation,rotation_mode="anchor"
            ha=x1_alignment, va="center",rotation_mode="anchor"
        )
        
        offset1 += size + PAD

    range_max = int(round(np.max(VALUES)/10)*10)
    HANGLES = np.linspace(0, 2 * np.pi, 200)
    for c in range(10,range_max+5,10):
        ax.plot(HANGLES, np.repeat(c + PLUS, 200), color= "#bebebe", lw=0.7)
    ylabel = np.arange(0.1, range_max/100+0.05, 0.1)
    ylabel = np.round(ylabel,2)
    ax.set_yticks(range(10+PLUS,range_max+PLUS+5,10))
    ax.set_yticklabels(list(map(str, ylabel)), color="grey", size=12)
    ax.set_rlabel_position(5)
    
    legend_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10) for color in palette]
    # add index to legend
    group_names_full = [str(i+1) + ". " + group_names_full[i] for i in range(len(group_names_full))]
    ax.legend(legend_handles,group_names_full,
        frameon=False,
        loc='upper center', 
        bbox_to_anchor=(0.5, 0), 
        ncol=3,
        fontsize=16,
        handlelength=1,
        handleheight=1,
        handletextpad=0.2,columnspacing=0.7)

    return fig
    
def circular_subplot(ax,VALUES, LABELS,PVALUE,COLORS,Group,p_threshold=0.05):
    # All this part is like the code above
    #VALUES = VALUES * 5
    #sort VALUES from low to high and sort LABELS accordingly
    idx = np.argsort(VALUES)
    VALUES = VALUES[idx]
    LABELS = LABELS[idx]
    PVALUE = PVALUE[idx]

    LABELS[PVALUE<=p_threshold] = LABELS[PVALUE<=p_threshold] + "*"

    ANGLES_N = len(VALUES)
    ANGLES = np.linspace(0, 2 * np.pi, num=ANGLES_N, endpoint=False)
    #WIDTH = (2 * np.pi) / len(ANGLES)

    GROUPS_SIZE = len(VALUES)

    PLUS = 15
    offset = 0
    OFFSET = np.pi / 2
    
    ax.set_theta_offset(OFFSET)
    ax.set_ylim(0, np.max(VALUES) + PLUS + 10)
    ax.set_frame_on(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    HANGLES = np.linspace(0, 2 * np.pi, 200)
    ax.plot(HANGLES, np.repeat(PLUS, 200), color= COLORS, lw=0.7)

    range_max = int(round(np.max(VALUES)/10)*10)
    ax.vlines(
        ANGLES, PLUS,range_max +  PLUS, color="#bebebe", linewidth=0.7
        #edgecolor="white", linewidth=2
    )
    for c in range(10,range_max+5,10):
        ax.plot(HANGLES, np.repeat(c + PLUS, 200), color= "#bebebe", lw=0.7)
    ylabel = np.arange(0.1, range_max/100+0.05, 0.1)
    ylabel = np.round(ylabel,2)
    ax.set_yticks(range(10+PLUS,range_max+PLUS+5,10))
    ax.set_yticklabels(list(map(str, ylabel)), color="grey", size=8)
    ax.fill(ANGLES, VALUES+PLUS, color=COLORS, alpha=0.4)
    add_labels_equalsize(ANGLES, np.max(VALUES) +PLUS, LABELS, OFFSET-np.pi / 20, ax,18)

    ax.vlines(ANGLES, PLUS,VALUES +  PLUS, color=COLORS, linewidth=1)
    #ax.plot(ANGLES, VALUES+PLUS, color=COLORS, linewidth=2, linestyle='solid', alpha=0.8)
    ax.scatter(ANGLES, VALUES +  PLUS , s=scale_to_interval(VALUES),color=COLORS)
    #add_labels_equalsize(ANGLES, np.max(VALUES) + PLUS, LABELS, OFFSET, ax,10)
    ax.text(
        x=0.5, y=0.5, s=Group, color="black", va="center", ha="center",
        ma="center", fontsize=16, fontweight="bold",
        linespacing=0.87, transform=ax.transAxes
    )
    #ax.set_title(Group, size=11, color="black",y=0.4)
    #plt.tight_layout()

def circular_singleplot(df,Group,group_names=None):
  # All this part is like the code above
    df_curr = df[df["group"]==Group]
    VALUES = df_curr["dice"].values * 100
    LABELS = df_curr["name"].values
    PVALUE = df_curr["p_value"].values

    fig, ax = plt.subplots(figsize=(3, 3), subplot_kw={"projection": "polar"})
    fig.patch.set_facecolor('white')
    if group_names is None:
        COLORS = "#1f77b4"
    else:
        palette = sns.color_palette(cc.glasbey, n_colors=len(group_names))
        # find index of group
        group_names = list(group_names)
        idx = group_names.index(Group)
        # COLORS is the i-th color in palette, i is the index of the group
        COLORS = palette[idx]
    
    circular_subplot(ax,VALUES, LABELS,PVALUE,COLORS,1)

    
def circular_atlases_all(df,group_names,p_threshold=0.05):
    if len(group_names) == 2:
        sub_row = 1
        sub_col = 2
    elif len(group_names) == 3:
        sub_row = 1
        sub_col = 3
    elif len(group_names) == 4:
        sub_row = 2
        sub_col = 2
    elif len(group_names) <= 6:
        sub_row = 2
        sub_col = 3
    elif len(group_names) <= 8:
        sub_row = 2
        sub_col = 4
    elif len(group_names) == 9:
        sub_row = 3
        sub_col = 3
    elif len(group_names) <= 12:
        sub_row = 3
        sub_col = 4
    elif len(group_names) <= 16:
        sub_row = 4
        sub_col = 4
    elif len(group_names) == 20:
        sub_row = 5
        sub_col = 4
    elif len(group_names) <= 24:
        sub_row = 6
        sub_col = 4

    palette = sns.color_palette(cc.glasbey, n_colors=len(group_names))
    # find index of group
    group_names = list(group_names)

    fig, axes = plt.subplots(sub_row, sub_col, figsize=(5 * sub_col, 5 * sub_row),subplot_kw={"projection": "polar"})
    fig.patch.set_facecolor('white')

    axes = axes.flatten()
    for Group in group_names:
        idx = group_names.index(Group)
        # COLORS is the i-th color in palette, i is the index of the group
        COLORS = palette[idx]
        
        df_curr = df[df["group"]==Group]
        VALUES = df_curr["dice"].values * 100
        LABELS = df_curr["name"].values
        PVALUE = df_curr["p_value"].values   

        ax = axes[idx]

        circular_subplot(ax,VALUES, LABELS,PVALUE,COLORS,idx+1)
    fig.subplots_adjust(wspace=2, hspace=1.1)
    return fig

def table_summary_plot(df,atlas_names_list, atlas_names_list_full):
    fig = {}
    c=0
    for i in range(0, len(atlas_names_list), 3):
        atlas_names = atlas_names_list[i:i+3]
        atlas_names_fullname = atlas_names_list_full[i:i+3]
        fig[c] = table_summary_singleplot(df,atlas_names,atlas_names_fullname)
        c+=1
    return fig

def table_summary_singleplot(df,atlas_names,atlas_names_fullname):
    for atlas_name in atlas_names:
        df_curr = df[df["group"]==atlas_name]
        df_curr = df_curr[['name','dice','p_value']]
        df_curr = df_curr.rename(columns={"dice": atlas_name + "_dice"})
        df_curr = df_curr.rename(columns={"name": atlas_name + "_name"})
        df_curr = df_curr.rename(columns={"p_value": atlas_name + "_p_value"})
        df_curr = df_curr.reset_index(drop=True)

        if atlas_name == atlas_names[0]:
            df_all = df_curr
        else:
            df_all = pd.concat([df_all,df_curr],axis=1)

    fig, ax = plt.subplots(figsize=(14, 8))
    #dice_list = ["UKBICA","MG360J12"]
    #name_list = ["name1","name2"]
    dice_list = [atlas_name + "_dice" for atlas_name in atlas_names]
    name_list = [atlas_name + "_name" for atlas_name in atlas_names]
    p_value_list = [atlas_name + "_p_value" for atlas_name in atlas_names]
    for name_c in name_list:
        df_all[name_c] = df_all[name_c].fillna('-')
    for dice_c in dice_list:
        # replace 0 with 0.0001
        df_all[dice_c] = df_all[dice_c].replace(0,0.0001)
        df_all[dice_c] = df_all[dice_c].fillna(0)
    for p_value_c in p_value_list:
        df_all[p_value_c] = df_all[p_value_c].fillna(0)

    col_defs = (
        [ColumnDefinition(
            name=name_col,
            title="Name",
            #width=1,
            textprops={"ha": "center"},
            border="left",
            group=atlas_names_fullname[name_list.index(name_col)],
        )
        for name_col in name_list]
        +
        [ColumnDefinition(
            name=dice_col,
            title="Dice",
            width=0.5,
            #plot_fn=bar,
            #plot_kw={
                    #"is_pct": True,
            #        "cmap": matplotlib.cm.magma,
            #        "plot_bg_bar": False,
            #        "annotate": True,
            #        "height": 0.5,
            #        "lw": 0.5,
            #        "formatter":dice_percent
            #        },
            #formatter=lambda x: round(x, 1)),
            formatter=dice_format,
            #textprops={
            #    "ha": "center",
            #    "bbox": {"boxstyle": "circle", "pad": 0.1},
            #},
            #cmap=normed_cmap(df_all[dice_col], cmap=matplotlib.cm.viridis, num_stds=4),
            group=atlas_names_fullname[dice_list.index(dice_col)],
        )
        for dice_col in dice_list]
        +
        [ColumnDefinition(
            name=p_value_col,
            title="P value",
            width=0.7,
            formatter=pvalue_format,
            textprops={
                "ha": "center",
                #"bbox": {"boxstyle": "circle", "pad": 0.1},
            },
            border="right",
            #cmap=normed_cmap(df_all[p_value_col], cmap=matplotlib.cm.YlGnBu, num_stds=4),
            group=atlas_names_fullname[p_value_list.index(p_value_col)],
        )
        for p_value_col in p_value_list]

    )
    # loop through every three columns
    df_all = df_all.set_index(name_list[0])
    tab = Table(df_all, 
                textprops={"ha": "center", "fontsize": 12},
                column_definitions=col_defs,
                row_dividers=True,
                footer_divider=True,
                row_divider_kw={"linewidth": 1, "linestyle": (0, (1, 5))},
                col_label_divider_kw={"linewidth": 1, "linestyle": "-"},
                column_border_kw={"linewidth": 1, "linestyle": "-"},
    ).autoset_fontcolors(colnames=dice_list + p_value_list)
    
    return fig

def draw_pvalue_mat(overlap_data, ref_atlas_name, other_atlas_name, minv, maxv, figfile):
    """
    Draw a given overlap matrix. 
    """
    overlap_data_reorder = copy.deepcopy(overlap_data)
    order_idx = pair_match(overlap_data_reorder)
    overlap_data = overlap_data[:, order_idx]
    reference_name = construct_network_name(ref_atlas_name)
    other_name = construct_network_name(other_atlas_name)

    other_name_new = [other_name[i] for i in order_idx]

    _, axes = plt.subplots(figsize=(8,8))
    sns.set(font_scale=1)
    sns.set_style("ticks")
    im = sns.heatmap(
            ax=axes,
            data=overlap_data,
            cmap="rocket",
            square=True,
            vmin=minv,
            vmax=maxv,
            linewidth=0.05,
            xticklabels=other_name_new,
            yticklabels=reference_name,
            cbar=False
            )

    frame_len=3.5
    axes.axhline(y=0, color='k',linewidth=frame_len)
    axes.axhline(y=overlap_data.shape[0], color='k',linewidth=frame_len)
    axes.axvline(x=0, color='k',linewidth=frame_len)
    axes.axvline(x=overlap_data.shape[1], color='k',linewidth=frame_len)

    axes.xaxis.tick_top()
    axes.tick_params(axis='x', which='major',labelsize=14,labelrotation=90,pad=10)
    axes.tick_params(axis='y', which='major',labelsize=14,pad=10)
        
    mappable = im.get_children()[0]
    cb=plt.colorbar(mappable,
                ax = axes,
                orientation = 'horizontal',
                location = 'bottom',
                label='Dice overlap',
                aspect=35.5,
                fraction=0.025,
                pad=0.08
                )
    cb.outline.set_linewidth(1.5)
    cb.outline.set_color('black')
    cb.ax.locator_params(nbins=10)
    cb.ax.tick_params(labelsize=10)
    cb.ax.tick_params(width=1.5)
    
    plt.tight_layout()
    plt.ioff()
    plt.savefig(figfile, dpi=200, bbox_inches='tight')

def table_summary_atlas(df):

    formatted_df = df
    other_name = df.columns.tolist()

    # Create a plot with a table
    fig, ax = plt.subplots(figsize=(24, 14))

    # Hide the axes
    ax.axis('off')

    # Plot the table with rotated column headers
    tab = table(ax, formatted_df, loc='center', colWidths=[0.1] * len(formatted_df.columns), cellLoc='center', rowLoc='center', cellColours=None)

    for key, cell in tab.get_celld().items():
        cell.set_height(0.06)
        cell.get_text().set_fontsize(16)
        cell.set_linewidth(2) 
    
    for i, key in enumerate(tab.get_celld().keys()):
        cell = tab[key]
        if key[0] == 0:
            cell.set_text_props(weight='bold', rotation=0, va='center', ha='center')
        
        if key[1] == -1:  # Check if it's a row name cell
            cell.set_text_props(weight='bold')
    plt.show()
    return fig

def report_data_network_correspondence(df,out_dir,atlas_names_list=None, atlas_names_list_full=None):
    # if directory doesn't exist, create it
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if atlas_names_list is None:
        atlas_names_list = df['group'].unique()
    if atlas_names_list_full is None:
        atlas_names_list_full = atlas_names_list            
    
    fig_chart = circular_barchart(df,atlas_names_list,atlas_names_list_full)
    #save figure
    #plt.savefig(fig_dir + "/circular_barchart.png", dpi=300, bbox_inches='tight')
    fig_spider = circular_atlases_all(df,atlas_names_list)
    #save figure
    #plt.savefig(fig_dir + "/circular_atlases_all.png", dpi=300, bbox_inches='tight')
    figs_table = table_summary_plot(df,atlas_names_list,atlas_names_list_full) 
    fig_chart.savefig(out_dir + "/circular_barchart.png", dpi=300, bbox_inches='tight')
    fig_spider.savefig(out_dir + "/circular_atlases_all.png", dpi=300, bbox_inches='tight')
    for key in figs_table:
        figs_table[key].savefig(out_dir + "/table_summary_plot" + str(key) + ".png", dpi=300, bbox_inches='tight')
    # export to csv
    df.to_csv(out_dir + "/network_correspondence.csv")

def report_atlas_network_correspondence(combined_data,ref_atlas_name,out_dir):
    # if directory doesn't exist, create it
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    reference_name = vis.construct_network_name(ref_atlas_name)
    iter_data = iter(combined_data.items())
    firstkey , firstvalue = next(iter_data)
    data_dim = firstvalue['Dice'].shape[0]
    # if reference_name length is 1, add index to reference_name
    if len(reference_name) == 1:
        reference_name = [reference_name[0] + str(i+1) for i in range(data_dim)]
    vis.draw_overlap_mat_all(combined_data, reference_name, figfile=out_dir)
    #for each key in combined_data, generate a table_summary_plot
    for key in combined_data:
        curr_data = combined_data[key]
        curr_dice = curr_data['Dice']
        curr_pval = curr_data['P value']
        
        other_name = vis.construct_network_name(key)
        new_array = np.array([f'{a:.4f}\n(p={b:.4f})' for a, b in zip(curr_dice.flatten(), curr_pval.flatten())]).reshape(curr_dice.shape)
        df = pd.DataFrame(new_array, index=reference_name, columns=other_name)
        fig_table = table_summary_atlas(df)
        
        fig_table.savefig(out_dir + "/table_summary_plot_" + key + ".png", dpi=300, bbox_inches='tight')
        # export to csv
        df.to_csv(out_dir + "/network_correspondence_" + key + ".csv")
