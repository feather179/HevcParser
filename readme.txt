基本是对照着 HEVC Spec 文档来写的逻辑，没有做特别的优化
实现了基础的 HEVC 码流解析功能，还有一些特性没有实现
部分代码参考了现有的开源软件
    1. CabacReader 中管理所有 Context 的代码结构参考 ffmpeg
    2. DCT/DST 变换部分从官方的 HM 参考软件中拷贝过来的

编译:
    使用了 xmake 工具来管理工程，本地使用 MSVC 编译器测试
    配置: xmake f -m debug  或者 xmake f -m release 
    编译: xmake
