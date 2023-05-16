# python
from matplotlib import pyplot as plt

def make_spider(df,row,type_num,title,color):
    # Create a color palette:
    my_palette = plt.cm.get_cmap("Set2", len(df.index))

    # number of variable
    categories=list(df)[:-1]
    # print(categories)
    N = len(categories)

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Initialise the spider plot
    ax = plt.subplot(1,type_num+1,row+1, polar=True, )

    # If you want the first axis to be on top:
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    # Draw one axe per variable + add labels labels yet
    plt.xticks(angles[:-1], categories, color='black', size=10)
    # ax.set_xticklabels(rotation=45)

    # Draw ylabels
    ax.set_rlabel_position(0)
    plt.yticks([25,50,75], ["25","50","75"], color="grey", size=9)
    # plt.ylim(0,100)

    # Ind1
    values=df.loc[row].drop('B-type').values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, color=color, linewidth=2, linestyle='solid')
    ax.fill(angles, values, color=color, alpha=0.4)
    # Add a title
    plt.title(title, size=16, color='black', y=1.1)
