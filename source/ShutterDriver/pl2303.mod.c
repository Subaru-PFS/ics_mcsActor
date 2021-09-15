#include <linux/module.h>
#include <linux/vermagic.h>
#include <linux/compiler.h>

MODULE_INFO(vermagic, VERMAGIC_STRING);

__visible struct module __this_module
__attribute__((section(".gnu.linkonce.this_module"))) = {
	.name = KBUILD_MODNAME,
	.init = init_module,
#ifdef CONFIG_MODULE_UNLOAD
	.exit = cleanup_module,
#endif
	.arch = MODULE_ARCH_INIT,
};

static const struct modversion_info ____versions[]
__used
__attribute__((section("__versions"))) = {
	{ 0xb6fd8534, __VMLINUX_SYMBOL_STR(module_layout) },
	{ 0x3f787b15, __VMLINUX_SYMBOL_STR(usb_serial_deregister_drivers) },
	{ 0x26bf2ee6, __VMLINUX_SYMBOL_STR(usb_serial_register_drivers) },
	{ 0x677000fc, __VMLINUX_SYMBOL_STR(tty_insert_flip_string_fixed_flag) },
	{ 0x980e92c1, __VMLINUX_SYMBOL_STR(usb_serial_handle_sysrq_char) },
	{ 0x2def689b, __VMLINUX_SYMBOL_STR(__tty_insert_flip_char) },
	{ 0x445e65c9, __VMLINUX_SYMBOL_STR(tty_flip_buffer_push) },
	{ 0x124edb5e, __VMLINUX_SYMBOL_STR(usb_serial_generic_open) },
	{ 0xa5347d7b, __VMLINUX_SYMBOL_STR(usb_clear_halt) },
	{ 0xe9bcba71, __VMLINUX_SYMBOL_STR(kmem_cache_alloc_trace) },
	{ 0x54fe1637, __VMLINUX_SYMBOL_STR(kmalloc_caches) },
	{ 0xf2997713, __VMLINUX_SYMBOL_STR(tty_termios_hw_change) },
	{ 0x5dc92376, __VMLINUX_SYMBOL_STR(usb_kill_urb) },
	{ 0xa146c04d, __VMLINUX_SYMBOL_STR(usb_serial_generic_close) },
	{ 0x4f6b400b, __VMLINUX_SYMBOL_STR(_copy_from_user) },
	{ 0x4f8b5ddb, __VMLINUX_SYMBOL_STR(_copy_to_user) },
	{ 0xdb7305a1, __VMLINUX_SYMBOL_STR(__stack_chk_fail) },
	{ 0x9a753209, __VMLINUX_SYMBOL_STR(tty_encode_baud_rate) },
	{ 0x409873e3, __VMLINUX_SYMBOL_STR(tty_termios_baud_rate) },
	{ 0x37a0cba, __VMLINUX_SYMBOL_STR(kfree) },
	{ 0xd90d4d66, __VMLINUX_SYMBOL_STR(usb_control_msg) },
	{ 0x29761f99, __VMLINUX_SYMBOL_STR(usb_serial_handle_dcd_change) },
	{ 0x79bf9e04, __VMLINUX_SYMBOL_STR(usb_serial_handle_break) },
	{ 0x879db016, __VMLINUX_SYMBOL_STR(tty_kref_put) },
	{ 0xf8a8a5ed, __VMLINUX_SYMBOL_STR(tty_port_tty_get) },
	{ 0xa6bbd805, __VMLINUX_SYMBOL_STR(__wake_up) },
	{ 0x77c1c1bd, __VMLINUX_SYMBOL_STR(dev_err) },
	{ 0x65bc16ec, __VMLINUX_SYMBOL_STR(usb_submit_urb) },
	{ 0xcf9e89d1, __VMLINUX_SYMBOL_STR(__dynamic_dev_dbg) },
	{ 0x1916e38c, __VMLINUX_SYMBOL_STR(_raw_spin_unlock_irqrestore) },
	{ 0x680ec266, __VMLINUX_SYMBOL_STR(_raw_spin_lock_irqsave) },
	{ 0xb47d047e, __VMLINUX_SYMBOL_STR(current_task) },
	{ 0xbdfb6dbb, __VMLINUX_SYMBOL_STR(__fentry__) },
};

static const char __module_depends[]
__used
__attribute__((section(".modinfo"))) =
"depends=usbserial";

MODULE_ALIAS("usb:v067Bp2303d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v067Bp2304d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v067Bp04BBd*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v067Bp1234d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v067BpAAA2d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v067BpAAA0d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v067Bp0611d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v067Bp0612d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v067Bp0609d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v067Bp331Ad*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v067Bp0307d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v04BBp0A03d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v04BBp0A0Ed*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0557p2008d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0547p2008d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v056Ep5003d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v056Ep5004d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0EBAp1080d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0EBAp2080d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0DF7p0620d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0584pB000d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v2478p2008d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v1453p4026d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0731p0528d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v6189p2068d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v11F7p02DFd*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v04E8p8001d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v11F5p0001d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v11F5p0003d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v11F5p0004d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v11F5p0005d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v04A5p4027d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0745p0001d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v078Bp1234d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v10B5pAC70d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v079Bp0027d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0413p2101d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0E55p110Bd*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0731p2003d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v050Dp0257d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v058Fp9720d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v11F6p2001d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v07AAp002Ad*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v05ADp0FBAd*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v5372p2303d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v03F0p3524d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v04B8p0521d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v04B8p0522d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v054Cp0437d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v11ADp0001d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0B63p6530d*dc*dsc*dp*ic*isc*ip*in*");
MODULE_ALIAS("usb:v0B8Cp2303d*dc*dsc*dp*ic*isc*ip*in*");

MODULE_INFO(srcversion, "E164B23C8E8ECC4D7FCE30E");
