import threading


def worker(num):
    """这是线程执行的函数"""
    print(f"Worker {num}")

    # 返回值
    return num * 2


# 创建一个线程，并传递参数
t = threading.Thread(target=worker, args=(42,))

# 启动线程
t.start()

# 等待线程结束，并获取返回值
result = t.join()

print(f"Result: {result}")
