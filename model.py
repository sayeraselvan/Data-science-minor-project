import random
import numpy as np
import matplotlib.pyplot as plt
# complete script
# modify the below parameters to build a bacterial genomic model to construst a model with the below described parameters

repeats_percent = 0.5 #direct/inverted repeat percent
repeats = 2500 #number of repeat pairs
inter_replicore_on = 0.67 #among the repeats, the percentage of inter-replicore repeats
intra_replicore_on = 1 - inter_replicore_on # 1 - inter replicore repeats is the intra replicore repeats
genes = 4000 # number of genes
replicore = 0.8 #inital gene strand bias
generations = 10000 #number of evolutionary generations - given that a random repeat pair is chosen for inversion in every generation

Run = 17 #random number generator

my_str = []
repeats_str = []
repeats_str1 = []
repeats_str2 = []
gene_str = []
gene_str_1half = []
gene_str_2half = []
index_of_repeat = []
index_of_repeat1 = []
inter_replicore = []
intra_replicore = []
total_on_repeats = []
inter_on_repeats = []
intra_on_repeats = []
leading = []
lagging = []
first_element = 1
mid_element = int((2*repeats + genes)/2 + 1)
last_element = 2*repeats + genes

def make_str(type):
    for i in range(1, int(repeats)+1):
        my_test_str = 'A'+str(i)
        type.append(my_test_str)
    return type

make_str(repeats_str1)
make_str(repeats_str2)

num_to_swap = (1-repeats_percent) * repeats

for i in range(int(num_to_swap)):
    repeats_str2[i] = repeats_str2[i].swapcase()
    
def replicore1(type, type1, genes, replicore):
    for i in range(1, int(genes//2)+1):
        type.append(i)
    for i in range(0, int(genes*replicore)//2):
        type[i] = -abs(type[i]) 
    for i in range(int(genes//2)+1, int(genes)+1):
        type1.append(i)
    for i in range(0, int(genes*replicore)//2):
        type1[i] = -abs(type1[i])
    return type, type1

replicore1(gene_str_1half, gene_str_2half, genes, replicore)
random.shuffle(gene_str_1half)
random.shuffle(gene_str_2half)
# print(gene_str_1half, gene_str_2half)
# print(repeats_str1, repeats_str2)

sp1 = repeats - (repeats*repeats_percent*intra_replicore_on*0.5)
ep1 = repeats
sp2 = 0
ep2 = (repeats*repeats_percent*intra_replicore_on*0.5)
temp = repeats_str2[int(sp1):int(ep1)]  
repeats_str2[int(sp1):int(ep1)] = repeats_str1[int(sp2):int(ep2)]  
repeats_str1[int(sp2):int(ep2)] = temp 

sp3 = repeats - (repeats*repeats_percent*intra_replicore_on)
ep3 = repeats - (repeats*repeats_percent*intra_replicore_on*0.5)
sp4 = 0
ep4 = (repeats*repeats_percent*intra_replicore_on*0.5)
temp1 = repeats_str1[int(sp3):int(ep3)]  
repeats_str1[int(sp3):int(ep3)] = repeats_str2[int(sp4):int(ep4)]  
repeats_str2[int(sp4):int(ep4)] = temp1 

# print(repeats_str1, repeats_str2)

repeats_str = repeats_str1 + repeats_str2
half_replicore1 =  gene_str_1half + repeats_str1
half_replicore2 =  gene_str_2half + repeats_str2
# print(half_replicore1, half_replicore2)
# print(repeats_str)

indices1 = list(range(1, len(half_replicore1)))
indices2 = list(range(1, len(half_replicore2)))
random.shuffle(indices1)
random.shuffle(indices2)

half_replicore_1 = [half_replicore1[0]] + [half_replicore1[i] for i in indices1]
half_replicore_2 = [half_replicore2[0]] + [half_replicore2[i] for i in indices2]

my_str = half_replicore_1 + half_replicore_2
# print(my_str)

integers = [x for x in my_str if isinstance(x, int)]
first_index = 0
first_value = integers[first_index]
median_index = len(integers) // 2 
median_value = integers[median_index]
# last_index = len(integers) - 1
# last_value = integers[last_index]



def only_reverse(my_str,index_of_repeat):
    my_str[index_of_repeat[0]-1:index_of_repeat[1]] = my_str[index_of_repeat[0]-1:index_of_repeat[1]][::-1]
    return my_str

def swapcase_g_r(my_str,index_of_repeat):
    my_str[index_of_repeat[0]-1:index_of_repeat[1]-1] = [str(i).swapcase() if isinstance(i, str) and my_str[index_of_repeat[0]-1:index_of_repeat[1]].count(i) == 1 else -i if isinstance(i, int) and my_str[index_of_repeat[0]-1:index_of_repeat[1]].count(-i) == 1 else i if isinstance(i, str) else -i for i in my_str[index_of_repeat[0]-1:index_of_repeat[1]-1]]
    return my_str

def swapcase_r(my_str,index_of_repeat):
    my_str[index_of_repeat[0]-1:index_of_repeat[1]-1] = [i.swapcase() if isinstance(i, str) and my_str[index_of_repeat[0]-1:index_of_repeat[1]].count(i) == 1 else i for i in my_str[index_of_repeat[0]-1:index_of_repeat[1]-1]]
    return my_str

def checkpoints(my_str_1, my_str, first_value, mid_value, index_of_repeat):
    if first_value in my_str_1 or mid_value in my_str_1:
        swapcase_r(my_str,index_of_repeat)
        only_reverse(my_str,index_of_repeat)
        return my_str
    else:
        swapcase_g_r(my_str,index_of_repeat)
        only_reverse(my_str,index_of_repeat)
        return my_str
            
def indices(my_str, first_repeat):
    for i in range(len(my_str)):
        if my_str[i] == first_repeat[0]:  
            index_of_repeat.append(i+1)
            # or my_str[i] == first_repeat[0].swapcase()        
    return index_of_repeat

def points(my_str,first_value,median_value):
    first_point = None
    index_of_first_value = None
    for index, item in enumerate(my_str):
        if isinstance(item, int) and item == first_value:
            index_of_first_value = index
            break
    first_point = index_of_first_value + 1

    mid_point = None
    index_of_median = None
    for index, item in enumerate(my_str):
        if isinstance(item, int) and item == median_value:
            index_of_median = index
            break
    mid_point = index_of_median + 1
    
    return first_point, mid_point

def distributiongraph():
    x = []
    y = []
    y_1 = []
    y_2 = []

    for i in range(len(my_str)):
        if isinstance(my_str[i], str):
            if my_str.count(my_str[i]) == 1:
                y_1.append(y_1[-1]+1 if len(y_1) > 0 else 0)
            else:
                y_1.append(y_1[-1] if len(y_1) > 0 else 0)
        else:
            y_1.append(y_1[-1] if len(y_1) > 0 else 0)
                
    for i in range(len(my_str)):
        if isinstance(my_str[i], str):
            if my_str.count(my_str[i]) == 2:
                y_2.append(y_2[-1]+1 if len(y_2) > 0 else 0)
            else:
                y_2.append(y_2[-1] if len(y_2) > 0 else 0)
        else:
            y_2.append(y_2[-1] if len(y_2) > 0 else 0)
                
    for i in range(len(my_str)):
        if isinstance(my_str[i], int):
            y.append(y[-1]+1 if len(y) > 0 else 0)
        else:
            y.append(y[-1] if len(y) > 0 else 0)
        
    for i in range(1,len(my_str)+1):
        x.append(i)

    # print(x,y)
    # print(y_1,y_2)


    plt.plot(x, y_1, c='red', label='Directed repeats')
    plt.plot(x, y_2, c='green', label='Inverted repeats')
    plt.plot(x, y, c='blue', label = "Genes")
    plt.title('Distribution graph')
    plt.xlabel('Gene position')
    plt.ylabel('Count')
    plt.xlabel('Position')
    plt.ylabel('Count')
    plt.legend()
    plt.show()


    integers = [x for x in my_str if isinstance(x, int)]
    first_index = 0
    first_value = integers[first_index]
    median_index = len(integers) // 2 
    median_value = integers[median_index]

    first_point, mid_point = points(my_str, first_value, median_value)
    x_1 = []
    y_inter = []
    y_intra1 = []
    y_intra2 = []
    e = my_str[first_point-1:mid_point-1]
    f = my_str[mid_point-1:] + my_str[:first_point-1]
    elements_inter_replicore = []
    elements_intra_replicore1 = []
    elements_intra_replicore2 = []
            
    for g, i in zip(e, range(len(e))):
        if isinstance(g, str) and g in f:
            y_inter.append(y_inter[-1]+1 if len(y_inter) > 0 else 0)
        else:
            y_inter.append(y_inter[-1] if len(y_inter) > 0 else 0)        

    for g, i in zip(e, range(len(e))):
        if isinstance(g, str) and e.count(g) > 1:
            y_intra1.append(y_intra1[-1]+1/2 if len(y_intra1) > 0 else 0)
        else:
            y_intra1.append(y_intra1[-1] if len(y_intra1) > 0 else 0)  
            
    for h, i in zip(f, range(len(f))):
        if isinstance(h, str) and f.count(h) > 1:
            y_intra2.append(y_intra2[-1]+1/2 if len(y_intra2) > 0 else 0)
        else:
            y_intra2.append(y_intra2[-1] if len(y_intra2) > 0 else 0)
            
    # print(my_str)
    # print(y_inter,y_intra1,y_intra2)
    # print(len(y_inter),len(y_intra1),len(y_intra2))

    for i in range(1,int(len(my_str)/2 +1)):
        x_1.append(i)
        
    plt.plot(x_1, y_inter, c='red', label='Inter inverted repeats')
    plt.plot(x_1, y_intra1, c='green', label='Intra inverted repeats from ori to ter')
    plt.plot(x_1, y_intra2, c='blue', label='Intra inverted repeats from ter to ori')
    plt.title('Distribution graph')
    plt.xlabel('Inverted position')
    plt.ylabel('Count')
    plt.xlabel('Position')
    plt.ylabel('Count')
    plt.legend()
    plt.show()

# distributiongraph()


integers = [x for x in my_str if isinstance(x, int)]
first_index = 0
first_value = integers[first_index]
median_index = len(integers) // 2 
median_value = integers[median_index]

def points(my_str,first_value,median_value):
    first_point = None
    index_of_first_value = None
    for index, item in enumerate(my_str):
        if isinstance(item, int) and item == first_value:
            index_of_first_value = index
            break
    first_point = index_of_first_value + 1

    mid_point = None
    index_of_median = None
    for index, item in enumerate(my_str):
        if isinstance(item, int) and item == median_value:
            index_of_median = index
            break
    mid_point = index_of_median + 1
    
    return first_point, mid_point

for i in range(1,generations+1):
    x = random.randint(0, int(2*repeats)-1)
    first_repeat_l = repeats_str[x]
    # print(first_repeat_l)
    # print(str(i) + "generation")
    index_of_repeat=[]
    first_repeat = first_repeat_l.split()
    indices(my_str, first_repeat)
    print(first_repeat_l, index_of_repeat)
    s,t = 0,0
    
    first_point,mid_point = points(my_str,first_value,median_value)
    print(f"first point:{first_value} is {first_point} and mid_point:{median_value} is {mid_point}")
    
    if len(index_of_repeat) == 2:
        S1 = index_of_repeat[1] - index_of_repeat[0]
        S2 = index_of_repeat[0] + (last_element - index_of_repeat[1])
        if S2>=S1: 
            my_str1 = []
            my_str_1 = []
            my_str1 = my_str.copy()
            my_str_1 = my_str[index_of_repeat[0]-1: index_of_repeat[1]]
            checkpoints(my_str_1, my_str, first_value , median_value, index_of_repeat)
            
            first_point,mid_point = points(my_str,first_value,median_value)
            
            ori_ter = my_str[first_point-1:mid_point-1]
            ori_ter_count = 0                

            for num in ori_ter:     
                if isinstance(num, int):  
                    ori_ter_count += 1 
                    
            ter_ori = my_str[mid_point-1:] + my_str[:first_point-1]
            ter_ori_count = 0
            
            for num in ter_ori: 
                if isinstance(num, int): 
                    ter_ori_count += 1 
            
            oriter = ori_ter_count / genes
            terori = ter_ori_count / genes
            s = round(abs(oriter-terori), 3)
            if s < 0.25:
                print(f"{str(i)} + generation - ON stage (index is 2, straight, follows with {s})")
                
                lagging_genes = 0
                leading_genes = 0
                for num in my_str:
                    try:
                        num = int(num)
                        if num > 0:
                            lagging_genes += 1
                        elif num < 0:
                            leading_genes += 1
                    except ValueError:
                        continue
                leading.append(leading_genes)
                lagging.append(lagging_genes)
                
                a,b,c,d,k,j = 0,0,0,0,0,0
                e = my_str[first_point-1:mid_point-1]
                f = my_str[mid_point-1:] + my_str[:first_point-1]
                elements_inter_replicore = []
                elements_intra_replicore1 = []
                elements_intra_replicore2 = []
                            
                for g in e:
                    if isinstance(g, str) and g in f:
                        a += 1
                        elements_inter_replicore.append(g)
                    elif isinstance(g, str) and e.count(g) > 1:
                        b += 1/2
                        elements_intra_replicore1.append(g)
                    else: 
                        continue
                    
                for h in f:
                    if isinstance(h, str) and f.count(h) > 1:
                        c += 1/2
                        elements_intra_replicore2.append(h)
                    else:
                        continue
                d = (a + b + c) / repeats
                k = a/repeats
                j = (b+c)/repeats
                # print(f"Inter {a} elements are: {elements_inter_replicore}")
                # print(f"Intra {b+c} elements are: {elements_intra_replicore1+elements_intra_replicore2}")
                inter_replicore.append(a)
                intra_replicore.append(int(b+c))
                total_on_repeats.append(d)
                inter_on_repeats.append(k)
                intra_on_repeats.append(j)
                # print(f"{(i)} + generation + {my_str}")
                continue
                
            else: 
                print(str(i) + " generation - OFF stage (indexs is 2, straight, dont)")
                my_str = my_str1
                # print(my_str)
                
                first_point,mid_point = points(my_str,first_value,median_value)

                lagging_genes = 0
                leading_genes = 0
                for num in my_str:
                    try:
                        num = int(num)
                        if num > 0:
                            lagging_genes += 1
                        elif num < 0:
                            leading_genes += 1
                    except ValueError:
                        continue
                leading.append(leading_genes)
                lagging.append(lagging_genes)
                
                a,b,c,d,k,j = 0,0,0,0,0,0
                e = my_str[first_point-1:mid_point-1]
                f = my_str[mid_point-1:] + my_str[:first_point-1]
                elements_inter_replicore = []
                elements_intra_replicore1 = []
                elements_intra_replicore2 = []
                            
                for g in e:
                    if isinstance(g, str) and g in f:
                        a += 1
                        elements_inter_replicore.append(g)
                    elif isinstance(g, str) and e.count(g) > 1:
                        b += 1/2
                        elements_intra_replicore1.append(g)
                    else: 
                        continue
                    
                for h in f:
                    if isinstance(h, str) and f.count(h) > 1:
                        c += 1/2
                        elements_intra_replicore2.append(h)
                    else:
                        continue
                d = (a+b+c)/ repeats
                k = a/repeats
                j = (b+c)/repeats

                # print(f"Inter {a} elements are: {elements_inter_replicore}")
                # print(f"Intra {b+c} elements are: {elements_intra_replicore1+elements_intra_replicore2}")
                inter_replicore.append(a)
                intra_replicore.append(int(b+c))
                total_on_repeats.append(d)
                inter_on_repeats.append(k)
                intra_on_repeats.append(j)
                # print(f"{(i)} + generation + {my_str}")
                continue
        
        else:
            my_str2 = []
            my_str2 = my_str.copy()
            first_point,mid_point = points(my_str,first_value,median_value)
            
            list1 = my_str[index_of_repeat[1]-1:last_element]
            list2 = my_str[0:index_of_repeat[0]]
        
            list_combined = list1 + list2
            index_of_repeat1 = [1,len(list_combined)]
            
            checkpoints(list_combined, list_combined, first_value, median_value, index_of_repeat1)
            
            split_index = last_element-index_of_repeat[1]+1

            first_half = list_combined[:split_index]
            second_half = list_combined[split_index:]
            my_str[index_of_repeat[1]-1:last_element] = first_half
            my_str[0:index_of_repeat[0]] = second_half
            
            first_point,mid_point = points(my_str,first_value,median_value)
            
            ori_ter1 = my_str[first_point-1:mid_point-1]
            ori_ter1_count = 0                

            for num in ori_ter1:     
                if isinstance(num, int):  
                    ori_ter1_count += 1 
                    
            ter_ori1 = my_str[mid_point-1:] + my_str[:first_point-1]
            ter_ori1_count = 0
            
            for num in ter_ori1: 
                if isinstance(num, int): 
                    ter_ori1_count += 1 
            
            oriter1 = ori_ter1_count / genes
            terori1 = ter_ori1_count / genes
            t = round(abs(oriter1-terori1), 3)
            if t < 0.25:
                print(f"{str(i)} +  generation - ON stage (index is 2, opposite, follows with {t})")
                
                lagging_genes = 0
                leading_genes = 0
                for num in my_str:
                    try:
                        num = int(num)
                        if num > 0:
                            lagging_genes += 1
                        elif num < 0:
                            leading_genes += 1
                    except ValueError:
                        continue
                leading.append(leading_genes)
                lagging.append(lagging_genes)
                
                a,b,c,d,k,j = 0,0,0,0,0,0
                e = my_str[first_point-1:mid_point-1]
                f = my_str[mid_point-1:] + my_str[:first_point-1]
                elements_inter_replicore = []
                elements_intra_replicore1 = []
                elements_intra_replicore2 = []
                            
                for g in e:
                    if isinstance(g, str) and g in f:
                        a += 1
                        elements_inter_replicore.append(g)
                    elif isinstance(g, str) and e.count(g) > 1:
                        b += 1/2
                        elements_intra_replicore1.append(g)
                    else: 
                        continue
                    
                for h in f:
                    if isinstance(h, str) and f.count(h) > 1:
                        c += 1/2
                        elements_intra_replicore2.append(h)
                    else:
                        continue
                d = (a+b+c)/ repeats
                k = a/repeats
                j = (b+c)/repeats

                # print(f"Inter {a} elements are: {elements_inter_replicore}")
                # print(f"Intra {b+c} elements are: {elements_intra_replicore1+elements_intra_replicore2}")
                inter_replicore.append(a)
                intra_replicore.append(int(b+c))
                total_on_repeats.append(d)
                inter_on_repeats.append(k)
                intra_on_repeats.append(j)
                # print(f"{(i)} + generation + {my_str}")
                continue
                
            else: 
                print(str(i) + " generation - OFF stage (index is 2, opposite, dont)")
                my_str = my_str2
                # print(my_str)
                
                first_point,mid_point = points(my_str,first_value,median_value)

                lagging_genes = 0
                leading_genes = 0
                for num in my_str:
                    try:
                        num = int(num)
                        if num > 0:
                            lagging_genes += 1
                        elif num < 0:
                            leading_genes += 1
                    except ValueError:
                        continue
                leading.append(leading_genes)
                lagging.append(lagging_genes)
                
                a,b,c,d,k,j = 0,0,0,0,0,0
                e = my_str[first_point-1:mid_point-1]
                f = my_str[mid_point-1:] + my_str[:first_point-1]
                elements_inter_replicore = []
                elements_intra_replicore1 = []
                elements_intra_replicore2 = []
                            
                for g in e:
                    if isinstance(g, str) and g in f:
                        a += 1
                        elements_inter_replicore.append(g)
                    elif isinstance(g, str) and e.count(g) > 1:
                        b += 1/2
                        elements_intra_replicore1.append(g)
                    else: 
                        continue
                    
                for h in f:
                    if isinstance(h, str) and f.count(h) > 1:
                        c += 1/2
                        elements_intra_replicore2.append(h)
                    else:
                        continue
                d = (a+b+c)/ repeats
                k = a/repeats
                j = (b+c)/repeats

                # print(f"Inter {a} elements are: {elements_inter_replicore}")
                # print(f"Intra {b+c} elements are: {elements_intra_replicore1+elements_intra_replicore2}")
                inter_replicore.append(a)
                intra_replicore.append(int(b+c))
                total_on_repeats.append(d)
                inter_on_repeats.append(k)
                intra_on_repeats.append(j)
                # print(f"{(i)} + generation + {my_str}")
                continue
    else:
        print(str(i) + " generation - OFF stage (index - 1 or 0)")
        # print(my_str)

        lagging_genes = 0
        leading_genes = 0
        for num in my_str:
            try:
                num = int(num)
                if num > 0:
                    lagging_genes += 1
                elif num < 0:
                    leading_genes += 1
            except ValueError:
                continue
        leading.append(leading_genes)
        lagging.append(lagging_genes)
        
        a,b,c,d,k,j = 0,0,0,0,0,0
        e = my_str[first_point-1:mid_point-1]
        f = my_str[mid_point-1:] + my_str[:first_point-1]
        elements_inter_replicore = []
        elements_intra_replicore1 = []
        elements_intra_replicore2 = []
        
        for g in e:
            if isinstance(g, str) and g in f:
                a += 1
                elements_inter_replicore.append(g)
            elif isinstance(g, str) and e.count(g) > 1:
                b += 1/2
                elements_intra_replicore1.append(g)
            else: 
                continue
            
        for h in f:
            if isinstance(h, str) and f.count(h) > 1:
                c += 1/2
                elements_intra_replicore2.append(h)
            else:
                continue
        d = (a+b+c)/ repeats
        k = a/repeats
        j = (b+c)/repeats

        # print(f"Inter {a} elements are: {elements_inter_replicore}")
        # print(f"Intra {b+c} elements are: {elements_intra_replicore1+elements_intra_replicore2}")
        inter_replicore.append(a)
        intra_replicore.append(int(b+c))
        total_on_repeats.append(d)
        inter_on_repeats.append(k)
        intra_on_repeats.append(j)
        # print(f"{(i)} + generation + {my_str}")
        continue
               
print(f"Leading gene are: {leading}")
print(f"Lagging gene are: {lagging}")
print(f"Inter replicore ON are: {inter_replicore}")
print(f"Intra replicore ON are: {intra_replicore}")
print(f"Total ON repeats (a+b+c)/repeats are: {total_on_repeats}")
print(f"Inter ON percent are: {inter_on_repeats}")
print(f"Intra ON percent are: {intra_on_repeats}")


def graph():
    x_values_1 = []
    for i in range(1, generations+1):
        x_values_1.append(i)

    plt.plot(x_values_1,leading, label='Leading strand genes')
    plt.plot(x_values_1,lagging, label='Lagging starnd genes')
    plt.plot(x_values_1,inter_replicore, label='ON Inter-replicore')
    plt.plot(x_values_1,intra_replicore, label='ON Intra-replicore')
    plt.xlabel('Generation')
    plt.ylabel('Genes & Repeats')
    plt.title(f'Plot btw Gen, Genes & Repeats {Run}')
    plt.legend()
    plt.show()


def graph_repeats():
    x_values_2 = []
    for i in range(1, generations+1):
        x_values_2.append(i)

    plt.plot(x_values_2,inter_replicore, label='ON Inter-replicore')
    plt.plot(x_values_2,intra_replicore, label='ON Intra-replicore')
    plt.xlabel('Generation')
    plt.ylabel('Inter/Intra ON repeats')
    plt.title(f'Plot btw Gen, Inter/Intra ON repeats {Run}')
    plt.legend()
    plt.show()
    
def graph_genes():
    x_values_3 = []
    for i in range(1, generations+1):
        x_values_3.append(i)

    plt.plot(x_values_3,leading, label='Leading strand genes')
    plt.plot(x_values_3,lagging, label='Lagging starnd genes')

    plt.xlabel('Generation')
    plt.ylabel('Strand Bias')
    plt.title(f'Plot btw Gen, Strand Bias {Run}')
    plt.legend()
    plt.show()

def graph_total_on_repeats():
    x_values_4 = []
    for i in range(1, generations+1):
        x_values_4.append(i)

    plt.plot(x_values_4,total_on_repeats, label='Inter + intra on percent')
    plt.plot(x_values_4,inter_on_repeats, label='Inter on percent')
    plt.plot(x_values_4,intra_on_repeats, label='Intra on percent')
    plt.xlabel('Generation')
    plt.ylabel('ON repeats percent')
    plt.title(f'Plot btw Gen, ON repeats percent {Run}')
    plt.ylim(0, 1)
    plt.legend()
    plt.show()


leading1 = []
lagging1 = []
leading1 = leading.copy()
lagging1 = lagging.copy()

for i in range(len(leading1)):
    leading1[i] = leading1[i] / genes
for i in range(len(lagging1)):
    lagging1[i] = lagging1[i] / genes

def graph_gene_repeats_percent():
    x_values_5 = []
    for i in range(1, generations+1):
        x_values_5.append(i)

    plt.plot(x_values_5,leading1, label='Leading strand gene bias')
    plt.plot(x_values_5,lagging1, label='Lagging strand gene bias')
    plt.plot(x_values_5,inter_on_repeats, label='ON Inter-replicore')
    plt.plot(x_values_5,intra_on_repeats, label='ON Intra-replicore')
    # plt.plot(x_values_5,total_on_repeats, label='Inter + intra on percent')
    plt.xlabel('Generation')
    plt.ylabel('Genes & Repeats')
    plt.title(f'Plot btw Gen, Genes & Repeats {Run}')
    plt.ylim(0, 1)
    plt.legend()
    plt.show()

def graph_genes_percent():
    x_values_6 = []
    for i in range(1, generations+1):
        x_values_6.append(i)

    plt.plot(x_values_6,leading1, label='Leading strand gene bias')
    plt.plot(x_values_6,lagging1, label='Lagging strand gene bias')
    plt.xlabel('Generation')
    plt.ylabel('Gene strand bias')
    plt.title(f'Plot btw Gen, Genes strand {Run}')
    plt.ylim(0, 1)
    plt.legend()
    plt.show()
     
# graph()
# graph_genes()
# graph_repeats()
graph_total_on_repeats()
graph_genes_percent()
graph_gene_repeats_percent()
        


    

