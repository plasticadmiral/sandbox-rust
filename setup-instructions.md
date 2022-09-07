# Installation
Use the following command to install rustup:

`curl https://sh.rustup.rs -sSf | sh`

# Directory structure
Rustup metadata and toolchains will be installed into the Rustup home directory, located at:

`/home/hari/.rustup`

This can be modified with the `RUSTUP_HOME` environment variable.


The Cargo home directory is located at:

`/home/hari/.cargo` 

This can be modified with the `CARGO_HOME` environment variable. 


The cargo, rustc, rustup and other commands will be added to Cargo's bin directory, located at: 

`/home/hari/.cargo/bin`


This path will then be added to your PATH environment variable by modifying the profile files located at:

`/home/hari/.profile`\
`/home/hari/.bashrc`\
`/home/hari/.zshenv`


You can uninstall at any time with `rustup self uninstall` and these changes will be reverted.
