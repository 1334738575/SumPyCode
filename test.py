import binascii

# 提供的十六进制字符串
hex_string = "49 40 15 f9 a3 5e 8b 22"

# 将十六进制字符串转换为字节数组
byte_array = bytearray.fromhex(hex_string)

print("Byte array:", byte_array)

# 如果知道具体的解码方式，可以在这里进行解码
# 例如，如果是 XOR 解码，需要提供 XOR 的密钥
# 这里仅仅是一个示例，实际解码需要根据 TightVNC 的实现来确定