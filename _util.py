

#%%
import numpy as np
from scipy import interpolate
import math
from scipy.linalg import eig
#%%

def print_tips(*obj):
    print(obj)

#%%

def eigensolve_close(H):
    '''
    get eigensolution of hamiltonian 'H'.
    '''
    vals, vecs = eig(H)    
    for i in range(len(vecs[:,1])):
        idx=np.append(range(i),(-abs(vecs[i,i:])).argsort()+i) if i>0 else (-abs(vecs[i,i:])).argsort()   
        vecs=vecs[:,idx]
        vals=vals[idx]
    return np.real(vals), vecs

def eigensolve_sort(H,ascending = True):
    '''
      get eigensolution of hamiltonian 'H', default ascending order is True.
      The return eigenenergies are in ascending order is ascending is True, else they will be is descending order.
    '''
    vals, vecs = eig(H)    
    if ascending:
        idx = vals.argsort()
    else:
        idx = vals.argsort()[::-1] 
    vals = vals[idx]
    vecs = vecs[:,idx]
    return np.real(vals), vecs


############################### deal with dictionary  ###############################################

def ValuesToDict(name,operator='aa',value=0):
    return {'name':name,'operator':operator,'value':value}

def get_IQ_dict(I,Q,LO_Freq):
    return {'I':I,'Q':Q,'lo_freq':LO_Freq}

def get_checked_dict(_dict):
    '''
        change all np.ndarray to list, because json file doesnt support np.ndarray
    '''
    for key,value in _dict.items():
        if isinstance(value,np.ndarray):
            _dict[key]=list(value)
        elif isinstance(value,dict):
            _dict[key]=get_checked_dict(value)
    return _dict

def get_dictionary_value(dict_,key_):
    try:
        value = dict_[key_]
    except KeyError:
        print(f'dictionary has no {key_}, set to None')
        value = None
    return value

def get_maximum_arraylen_in_dict(dict_,max_len=0):
    for value in dict_.values():
        if isinstance(value,(list,np.ndarray)):
            array_len = len(value)
            max_len = max(array_len,max_len)
        elif isinstance(value,dict):
            max_len = get_maximum_arraylen_in_dict(value,max_len)
    return max_len

def get_start_and_end_index_in_array(arr_value):
    arr_len = len(arr_value)
    static_value =  arr_value[0]
    start_index = 0
    end_index =1
    for i,value in enumerate(arr_value[1:]):
        if value != static_value:
            start_index = i
            break
    for i,value in enumerate(arr_value[::-1]):
        if value != static_value:
            end_index = arr_len - i + 1
            break
    return (max(0,start_index),min(end_index,arr_len))

def get_start_and_end_index_in_dict(dict_,start_index=None,end_index=None):
    for value in dict_.values():
        if isinstance(value,(list,np.ndarray)):
            arr_start_index,arr_end_index = get_start_and_end_index_in_array(value)
            if (not start_index) or (start_index > arr_start_index):
                start_index = arr_start_index
            if (not end_index) or (end_index < arr_end_index ):
                end_index = arr_end_index
        elif isinstance(value,dict):
            (start_index,end_index) = get_start_and_end_index_in_dict(value,start_index,end_index)
    return (start_index,end_index)

def cut_dict_by_start_and_end_index(dict_,start_index,end_index):
    if (start_index ==None) or (end_index ==None):
        return dict_
    for key,value in dict_.items():
        if isinstance(value,(list,np.ndarray)):
            try:
                dict_[key] = value[start_index:end_index]
            except Exception:
                dict_[key] = value[start_index:]
        elif isinstance(value,dict):
            dict_[key]=cut_dict_by_start_and_end_index(value,start_index,end_index)
    return dict_


#######################################################

def get_format_IQ(I,Q,lo_freq):
    I=np.asarray(I) if isinstance(I,(list,np.ndarray)) else np.asarray([I])
    Q=np.asarray(Q) if isinstance(Q,(list,np.ndarray)) else np.asarray([Q])
    if len(I)>len(Q):
        Q=np.append(Q,np.zeros(len(I)-len(Q)))
    elif len(I)<len(Q):
        I=np.append(I,np.zeros(len(Q)-len(I)))
    return I+1j*Q,lo_freq

def trunc_array(array,arr_len):
    ## if length of array > arr_len, trunc it
    ## else: repeat the last value

    if len(array)<arr_len:
        array=np.append(array,array[-1]*np.ones(arr_len-len(array)))
    elif len(array)>arr_len:
        array=array[0:arr_len]
    return array

def get_Interpolate_pulse(orginal_tlist,pulse,new_tlist=None,k=3):
    spline = interpolate.splrep(orginal_tlist,pulse,k=k)
    if new_tlist is None:
        return spline
    else:
        return interpolate.splev(new_tlist,spline)


################  str #################

def multiply_list(list_):
    product_ = 1
    for value in list_:
        product_ *= value
    return product_

def sstr2idx(sstr='000',dims_list=[3,3,4]):
    ## dims_list= [3,4,3]  ,  '010' -- 3
    inv_sstr = sstr[::-1]
    inv_dim = dims_list[::-1]
    idx = 0
    for i,value in enumerate(inv_sstr):
        idx += float(value) * multiply_list(inv_dim[0:i])
    return int(idx)

def idx2sstr(idx=0,dims_list=[3,3,4]):
    sstr=''
    mul = multiply_list(dims_list)
    for value in dims_list:
        mul /= value
        qs = idx // mul
        sstr += str(int(qs))
        idx -=  qs * mul
    return sstr

def split_str_func(func_str):
    str_copy=list(func_str)
    i=0
    for s in func_str:
        if (s=='+' or s=='-'):
            str_copy.insert(i,',')
            i+=1
        i+=1
    return ''.join(str_copy).split(',')

###################  level truncation #############

def trunc_to_lowest_two_level(matrix,nTrunc=3):
    if len(matrix)==1:
        matrix_type = 'bra'
        matrix_len = len(matrix[0])
    elif len(matrix[0])==1:
        matrix_type = 'ket'
        matrix_len =len(matrix)
    elif len(matrix) == len(matrix[0]):
        matrix_type ='oper'
        matrix_len = len(matrix)
    else:
        print('wrong dimension!')

    if isinstance(nTrunc,(int,float)):
        nTrunc = np.ones(int(matrix_len / nTrunc),dtype=int )

    del_idx=[]
    for i in range(matrix_len):
        sstr = idx2sstr(i,dims_list=nTrunc)
        if (sstr.count('0') + sstr.count('1')) < len(sstr):
            del_idx.append(i)

    for position_num in del_idx[::-1]:
        if matrix_type=='ket':
            matrix=np.delete(matrix,position_num,0)         
        elif matrix_type=='bra':
            matrix=np.delete(matrix,position_num,1)          
        elif matrix_type=='oper':
            matrix=np.delete(matrix,position_num,1)  
            matrix=np.delete(matrix,position_num,0)              
    return matrix

def trunc_to_specific_subspace(matrix,nTrunc=3,subspace=['000','001','100','101']):
    ## subspace
    if len(matrix)==1:
        matrix_type = 'bra'
        matrix_len = len(matrix[0])
    elif len(matrix[0])==1:
        matrix_type = 'ket'
        matrix_len =len(matrix)
    elif len(matrix) == len(matrix[0]):
        matrix_type ='oper'
        matrix_len = len(matrix)
    else:
        print('wrong dimension!')

    if isinstance(nTrunc,(int,float)):
        nTrunc = np.ones(int(matrix_len / nTrunc),dtype=int )

    qubit_num = len(idx2sstr(0,dims_list=nTrunc))
    if qubit_num != len(subspace[0]):
        raise Exception(f'Dimension of subspace:{subspace[0]} is not equal to qubits number: {qubit_num}' )

    del_idx=[]
    for i in range(matrix_len):
        sstr = idx2sstr(i,dims_list=nTrunc)
        if sstr not in subspace:
            del_idx.append(i)

    for position_num in del_idx[::-1]:
        if matrix_type=='ket':
            matrix=np.delete(matrix,position_num,0)         
        elif matrix_type=='bra':
            matrix=np.delete(matrix,position_num,1)          
        elif matrix_type=='oper':
            matrix=np.delete(matrix,position_num,1)  
            matrix=np.delete(matrix,position_num,0)              
    return matrix





#%%
if __name__ == '__main__':
    print(idx2sstr(4**6-1,[4,4,4,4,4,4]))
    print(sstr2idx('103',[3,3,4]))
