import math
import marcumq
import itertools

def calculate_var_psi(du, q):
    # 1. 计算核心参数 a
    denominator = 2 ** (du + 1)
    a = math.ceil(q / denominator)
    
    # 2. 中间值概率
    p_mid = (2 ** du) / q
    
    # 3. 计算中间值的平方和部分（1~a-1的平方和）
    if a - 1 >= 1:
        sum_x_sq = (a - 1) * a * (2 * a - 1) // 6
    else:
        sum_x_sq = 0  # a=1时，无中间值，平方和为0
    sum_sq_mid = 2 * sum_x_sq * p_mid
    
    # 4. 计算两端值的概率和平方和部分
    total_p_mid = (2 * a - 1) * p_mid
    p_end = (1 - total_p_mid) / 2
    sum_sq_end = 2 * (a ** 2) * p_end
    
    # 5. 总方差（E(X)=0，直接等于E(X²)）
    var_psi = sum_sq_mid + sum_sq_end
    return var_psi

def calculate_sigma_delta_squared(n,N,m,alpha1,alpha2,var_psi):
    # 分步计算（拆分公式，降低计算错误概率）
    term1 = (n * N * alpha1 ** 2) / 4          # 公式第一项
    term2 = (m * N * alpha1 / 2) * (alpha2 / 2 + var_psi)  # 公式第二项
    term3 = alpha2 / 2                         # 公式第三项
    sigma_result = term1 + term2 + term3       # 最终结果
    return math.sqrt(sigma_result)

def calculate_A(N,q,dv,sigma_delta):
    exponent = dv + 1                # 计算指数：dᵥ+1
    denominator_power = 2 ** exponent# 计算 2^(dᵥ+1)
    ceil_part = math.ceil(q / denominator_power)  # 向上取整 ⌈q/2^(dᵥ+1)⌉
    sqrt_N = math.sqrt(N)            # 计算 √N
    numerator = sqrt_N * ceil_part   # 计算分子：√N × 向上取整结果
    a = numerator / sigma_delta      # 最终计算 a 的值
    return a

def calculate_a(N, miu, sigma_delta):
    return math.sqrt(sum(x**2 for x in miu[:N]))/sigma_delta

def generate_all_miu(N, border):
    single_values = list(range(-border, border + 1))
    all_miu = [list(combination) for combination in itertools.product(single_values, repeat=N)]
    return all_miu

def calculate_b(Rp,tau_lambda2,sigma_delta):
    # 3. 核心公式计算
    numerator = (1 - 2**(-Rp)) * tau_lambda2  # 分子
    denominator = 2 * sigma_delta                  # 分母
    b = numerator / denominator
    return b

def calculate_Nu(q,dv):
    denominator = 2 ** (dv + 1)
    border = math.ceil(q / denominator)
    Nu = border*2+1
    return border, Nu


if __name__ == "__main__":
    N = 256
    q = 7681
    n = 3
    m = 3
    alpha1 = 3
    alpha2 = 2
    Rp = 2
    du = 11
    dv = 8

    tau = q/2
    var = calculate_var_psi(du, q)
    delta = calculate_sigma_delta_squared(n,N,m,alpha1,alpha2,var)
    A = calculate_A(N,q,dv,delta)
    B = calculate_b(Rp,tau,delta)
    error = marcumq.marcumq(nu=N/2,a=A,b=B)
    log_error = math.log2(error)
    print(log_error)
