const mongoose = require('mongoose');

const userSchema = new mongoose.Schema({
    name: String,
    orcid: String,
    role: {
        required: true,
        type: String,
        enum: ['visitor', 'curator'],
        default: 'visitor'
    }
  });
  
const UserModel = mongoose.model('users', userSchema);

module.exports = UserModel;