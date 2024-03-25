const mongoose = require('mongoose');

const userSchema = new mongoose.Schema({
    email: String,
    password: String,
    name: String,
    lastname: String,
    role: {
        required: true,
        type: String,
        enum: ['visitor', 'curator'],
        default: 'visitor'
    }
  });
  
const UserModel = mongoose.model('users', userSchema);

module.exports = UserModel;