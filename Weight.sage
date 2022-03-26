from collections import defaultdict 

class Weight:
    
    def __init__(self,root_system):
        """
        Initializes the object with a RootSystem
        example: input =  RootSystem(['A',5])
        """
        global R;
        
        R = root_system.root_lattice()
        self.__helper()
        
    def __helper(self):
        """
        computes h_order and alpha_sum
        """
        global ref,h_order,alpha_sum;
        
        W = WeylGroup(R,prefix="s")
        h_order = defaultdict(list) 
        alpha_sum = 0
        ref = W.reflections()
        for s in ref.keys():
            h_order[ref[s].length()].append(s)
            alpha_sum += s
       
            
    def calculate_weight(self,x):
        """
        calculates weight for an WeylGroup element
        """
        
        len_x = x.length()
        if len_x == 0:
            return 0
        for i in reversed(h_order):
            for alpha in h_order[i]:
                if ( ( x*ref[alpha]).length() == len_x-i ) and ( i == alpha_sum.scalar(alpha.associated_coroot()) - 1 ):
                    return alpha.associated_coroot()+self.weight(x*ref[alpha])
