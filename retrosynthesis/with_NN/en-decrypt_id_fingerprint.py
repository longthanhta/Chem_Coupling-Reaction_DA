## Encryption
def six_bits_to_code(x:str)->int:
    # size of x must be 6, 2^6 = 64 characters
    # note that the order of conversion is quite reversed
    return sum([int(x[i])*pow(2,i) for i in range(6)])

print(six_bits_to_code("000001"))
chars = "abcdefghijklnmopqrstuvwxyzABCDEFJHIJKLNMOPQRSTUVWXYZ0123456789$@"  
print(len(chars))

def encrypt_fp(fp: str)->str:
    
    """ Convert binary string to general string {0,1} -> {a:z,A:Z,0:9,$@} 
    (64 characters). Since we have 2048 bits, we ADD 4 0's at the end
    to take 342 6-bit phrases. When use CONVERT BACK to fp, PLEASE REMEMBER 
    to REMOVE 4 0's at the end"""
    
    # fp must have size 2048
    fp += "0000"
    out = ""
    for i in range(342):
        out += chars[six_bits_to_code(fp[i*6:i*6+6])]
    return out


## Decryption
def code_to_six_bit(n:int)->list:
    # 0 <= x <= 63
    return [(n//pow(2,i))%2 for i in range(6)]
reverse = dict()
for i in range(64):
    reverse[chars[i]] = code_to_six_bit(i)
def decrypt_fp(en_fp: str)->list: # input for NN model!
    ## Remember to remove 4 last 0's
    out = []
    for x in en_fp:
        out = out+reverse[x]
    return out[:-4]


fp = decrypt_fp('a'*341+'b')
print(len(fp))
print(fp)

''' 

## Conversion
file = open("id_fingerprint.txt",'r')
file.readline() # columns' title
out = open("encrypted_id_fp.csv",'w')

for line in file:
    try:
        [id_ ,fp_] = line[:-1].split(',')
        en_fp_ = encrypt_fp(fp_)
        out.write(id_+','+en_fp_+'\n')
    except:
        print(id_)

file.close()
out.close()

# It costs 1h38m to run
'''